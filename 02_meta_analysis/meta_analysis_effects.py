"""Dependence-aware exploratory meta-analysis of nutrient uptake lnRR effects.

The functions are imported by Meta_analysis_temperature_irradiance.ipynb.
Original workbook values are never modified.
"""

from pathlib import Path
import sys
from math import isfinite

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
from scipy.stats import t
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
from project_config import DATA_FILE, DATA_SHEET, RESULTS_DIR
SOURCE = DATA_FILE
OUTPUT = RESULTS_DIR / "meta_analysis"


def prepare(source=SOURCE):
    df = pd.read_excel(source, sheet_name=DATA_SHEET)
    required = ["Effect_ID", "Study_ID", "Experiment_ID", "Driver_harmonised",
                "Environmental_contrast_type", "Environmental_change_direction",
                "Nutrient_form", "Author_replication_effect_status",
                "Author_replication_lnRR", "Author_replication_sampling_variance_lnRR"]
    absent = sorted(set(required) - set(df.columns))
    if absent:
        raise ValueError(f"Missing columns: {absent}")
    if df.Effect_ID.isna().any() or df.Effect_ID.duplicated().any():
        raise ValueError("Effect_ID must be present and unique")
    eligible = df.Author_replication_effect_status.eq("Eligible")
    temp = df.Driver_harmonised.eq("Temperature") & df.Environmental_contrast_type.eq("Temperature level")
    light = df.Driver_harmonised.eq("Irradiance") & df.Environmental_contrast_type.isin(
        ["Irradiance intensity", "Light versus darkness"])
    pooling = df.Author_replication_pooling_status.fillna("").astype(str)
    selected = (pooling.str.startswith("Candidate after dependence-aware contrast selection") |
                pooling.str.startswith("Candidate after within-study aggregation") |
                pooling.str.startswith("Candidate within Vmax stratum") |
                pooling.str.startswith("Selected for primary driver-nutrient subset") |
                pooling.str.startswith("Eligible for the corresponding irradiance-by-nutrient pooling") |
                pooling.str.startswith("Different species and nutrients originate from one study"))
    out = df.loc[eligible & (temp | light) & selected].copy()
    out["y_original"] = pd.to_numeric(out.Author_replication_lnRR, errors="coerce")
    out["v"] = pd.to_numeric(out.Author_replication_sampling_variance_lnRR, errors="coerce")
    lower = out.Environmental_change_direction.isin(["Lower temperature", "Lower irradiance"])
    higher = out.Environmental_change_direction.isin(["Higher temperature", "Higher irradiance"])
    unknown = ~(lower | higher)
    if unknown.any():
        bad = out.loc[unknown, ["Effect_ID", "Environmental_change_direction"]]
        raise ValueError(f"Unorientable exposure direction:\n{bad.to_string(index=False)}")
    out["orientation_multiplier"] = np.where(lower, -1, 1)
    out["y"] = out.y_original * out.orientation_multiplier
    valid = np.isfinite(out.y) & np.isfinite(out.v) & (out.v > 0)
    if not valid.all():
        raise ValueError(f"Invalid lnRR or sampling variance for {out.loc[~valid, 'Effect_ID'].tolist()}")
    out["study_label"] = out.First_author.astype(str) + " (" + out.Year.astype(str) + ")"
    audit = pd.DataFrame({
        "group": ["all records", "author-eligible", "direct temperature/intensity/light-dark",
                  "eligible direct contrasts", "selected (no alternatives/unresolved/sensitivity)",
                  "analysis input"],
        "rows": [len(df), int(eligible.sum()), int((temp | light).sum()),
                 int((eligible & (temp | light)).sum()),
                 int((eligible & (temp | light) & selected).sum()), len(out)],
    })
    return out, audit


def _cov_block(group, rho, tau_study, tau_experiment):
    v = group.v.to_numpy(float)
    exp_id = group.Experiment_ID.astype(str).to_numpy()
    same_exp = (exp_id[:, None] == exp_id[None, :]).astype(float)
    sampling = np.diag(v) + rho * np.sqrt(np.outer(v, v)) * (same_exp - np.eye(len(v)))
    return sampling + tau_study * np.ones((len(v), len(v))) + tau_experiment * same_exp


def _gls_components(data, rho, tau_study, tau_experiment):
    a = b = logdet = quad = 0.0
    pieces = []
    for study_id, group in data.groupby("Study_ID", sort=False):
        y = group.y.to_numpy(float)
        cov = _cov_block(group, rho, tau_study, tau_experiment)
        try:
            cf = cho_factor(cov, lower=True, check_finite=False)
        except np.linalg.LinAlgError:
            return None
        one = np.ones(len(group))
        wy = cho_solve(cf, y, check_finite=False)
        w1 = cho_solve(cf, one, check_finite=False)
        ai = float(one @ w1)
        bi = float(one @ wy)
        a += ai
        b += bi
        logdet += 2 * np.log(np.diag(cf[0])).sum()
        quad += float(y @ wy)
        pieces.append((study_id, ai, bi))
    mu = b / a
    reml = logdet + np.log(a) + quad - b * b / a
    return mu, a, reml, pieces


def fit_meta(data, rho=0.0, start=None):
    """Two nested random intercepts; REML, with study-cluster robust mean SE.

    rho is an assumed *sampling-error* correlation within an experiment,
    not an estimated biological correlation. Run several rho values.
    """
    if not 0 <= rho < 1:
        raise ValueError("rho must lie in [0, 1)")
    if data.empty:
        return None
    k = data.Study_ID.nunique()
    n = len(data)
    # Variance components can legitimately be zero. Optimise on their original scale.
    scale = max(float(data.y.var()), float(data.v.median()), 0.01)
    def objective(x):
        result = _gls_components(data, rho, x[0], x[1])
        return 0.5 * result[2] if result is not None else 1e100
    starts = ([list(start)] if start is not None else
              [[0., 0.], [scale / 4, scale / 4], [scale, 0.], [0., scale]])
    fits = [minimize(objective, s, method="L-BFGS-B", bounds=[(0, None), (0, None)],
                     options={"maxiter": 300}) for s in starts]
    best = min(fits, key=lambda r: r.fun)
    if not np.isfinite(best.fun):
        raise RuntimeError("REML optimisation failed")
    tau_study, tau_exp = (max(float(z), 0) for z in best.x)
    mu, a, _, pieces = _gls_components(data, rho, tau_study, tau_exp)
    # Cluster sandwich reflects the fact that rows from one study are dependent.
    se_model = np.sqrt(1 / a)
    if k >= 2:
        score = sum((bi - ai * mu) ** 2 for _, ai, bi in pieces)
        se_cluster = np.sqrt(k / (k - 1) * score / a ** 2)
        se = max(se_cluster, se_model)
        critical = float(t.ppf(0.975, k - 1))
    else:
        se_cluster = np.nan
        se = np.nan
        critical = np.nan
    lo, hi = mu - critical * se, mu + critical * se
    return dict(effects=n, studies=k, experiments=data.Experiment_ID.nunique(),
                rho=rho, lnRR=mu, SE_model=se_model, SE_cluster=se_cluster,
                CI_low=lo, CI_high=hi, percent_change=100 * np.expm1(mu),
                percent_low=100 * np.expm1(lo) if isfinite(lo) else np.nan,
                percent_high=100 * np.expm1(hi) if isfinite(hi) else np.nan,
                tau2_study=tau_study, tau2_experiment=tau_exp,
                convergence=bool(best.success))


def run_models(data, rhos=(0, 0.5, 0.8), min_studies=3):
    rows = []
    for driver in ["Temperature", "Irradiance"]:
        subset = data.loc[data.Driver_harmonised.eq(driver)]
        groups = [("All nutrients (exploratory)", subset)] + list(subset.groupby("Nutrient_form"))
        for nutrient, group in groups:
            # Two-study models can be calculated but cannot support a stable synthesis.
            if group.Study_ID.nunique() < min_studies:
                continue
            for rho in rhos:
                fit = fit_meta(group, rho)
                rows.append(dict(Driver=driver, Nutrient=nutrient, **fit))
    return pd.DataFrame(rows)


def _display_points(data, rho, level):
    points = []
    key = "Experiment_ID" if level == "Experiment" else "Study_ID"
    for group_id, group in data.groupby(key, sort=False):
        cov = _cov_block(group, rho, 0., 0.)
        cf = cho_factor(cov, lower=True, check_finite=False)
        one = np.ones(len(group))
        a = float(one @ cho_solve(cf, one, check_finite=False))
        b = float(one @ cho_solve(cf, group.y.to_numpy(float), check_finite=False))
        if level == "Experiment":
            suffix = str(group_id).split("_")[-1]
            label = f"{group.study_label.iloc[0]} — {suffix}"
        else:
            label = group.study_label.iloc[0]
        points.append(dict(group_id=group_id, label=label,
                           effects=len(group), lnRR=b / a, SE=np.sqrt(1 / a)))
    return pd.DataFrame(points).sort_values(["label", "group_id"])


def forest_plot(data, model_row, output_dir=OUTPUT, display_level="Experiment"):
    """Plot experiment-level GLS summaries; pooled diamond uses all selected effects."""
    if display_level not in ("Experiment", "Study"):
        raise ValueError("display_level must be 'Experiment' or 'Study'")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    driver, nutrient, rho = model_row.Driver, model_row.Nutrient, model_row.rho
    group = data.loc[data.Driver_harmonised.eq(driver)]
    if nutrient != "All nutrients (exploratory)":
        group = group.loc[group.Nutrient_form.eq(nutrient)]
    points = _display_points(group, rho, display_level)
    count = len(points)
    fig, ax = plt.subplots(figsize=(9, max(4.5, 0.43 * (count + 3))))
    yy = np.arange(count, 0, -1)
    ci = 1.96 * points.SE.to_numpy()
    ax.errorbar(points.lnRR, yy, xerr=ci, fmt="o", color="#176B86",
                ecolor="#176B86", capsize=3, ms=5, linewidth=1)
    ax.axvline(0, color="#666666", linestyle="--", linewidth=1)
    ax.set_yticks(yy, [f"{r.label}  (m={r.effects})" for r in points.itertuples()])
    pooled_y = 0
    ax.plot(model_row.lnRR, pooled_y, "D", color="#9C3D10", ms=8)
    ax.plot([model_row.CI_low, model_row.CI_high], [pooled_y, pooled_y],
            color="#9C3D10", lw=2)
    ax.set_ylim(-1, count + 1)
    ax.set_xlabel("ln response ratio, oriented towards higher exposure")
    ax.set_title(f"{driver}: {nutrient} | points by {display_level.lower()} | assumed sampling ρ = {rho:g}")
    ax.set_yticks(list(yy) + [0],
                  [f"{r.label}  (m={r.effects})" for r in points.itertuples()] + ["REML pooled effect"])
    ax.grid(axis="x", color="#D9E1E4", linewidth=0.6)
    fig.tight_layout()
    filename = f"forest_{driver.lower()}_{str(nutrient).lower().replace(' ', '_').replace('(', '').replace(')', '')}_by_{display_level.lower()}_rho{str(rho).replace('.', 'p')}.png"
    path = output_dir / filename
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path
