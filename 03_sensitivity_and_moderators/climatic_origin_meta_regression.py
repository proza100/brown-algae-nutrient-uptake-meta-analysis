"""Exploratory climatic-origin meta-regressions for nitrate uptake.

Uses the same selected lnRR effects, nested study/experiment covariance,
REML variance components, and rho sensitivity assumptions as the primary
meta-analysis. Original workbook and analysis-input files are not modified.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
from scipy.stats import t


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
from project_config import RESULTS_DIR
INPUT = RESULTS_DIR / "meta_analysis" / "analysis_input_effects.csv"
OUTPUT = RESULTS_DIR / "meta_analysis"


COLD_STUDIES = {
    "Bommarito_2024", "Ding_2025", "Ding_2026", "Fales_2023",
    "Korb_2000", "Ozaki_2001", "Topinka_1978", "Umanzor_2023", "Xu_2011",
}
WARM_STUDIES = {
    "Almeida-Saa_2024", "Bonomi-Barufi_2025", "Gerard_1982", "Manley_1985",
    "SanchezBarredo_2020", "Sato_2023", "Smith_2021", "Umanzor_2021",
}


def add_climatic_origin(df):
    out = df.copy()
    out["Climatic_origin"] = np.select(
        [out.Citation_key.isin(COLD_STUDIES), out.Citation_key.isin(WARM_STUDIES)],
        ["Cold-origin", "Warm-origin"], default="Unclassified",
    )
    # Edwards contains two separately analysed populations.
    out.loc[out.Experiment_ID.eq("S008_E01"), "Climatic_origin"] = "Warm-origin"
    out.loc[out.Experiment_ID.eq("S008_E02"), "Climatic_origin"] = "Cold-origin"
    if out.Climatic_origin.eq("Unclassified").any():
        bad = sorted(out.loc[out.Climatic_origin.eq("Unclassified"), "Citation_key"].unique())
        raise ValueError(f"Unclassified studies: {bad}")
    out["warm"] = out.Climatic_origin.eq("Warm-origin").astype(float)
    return out


def covariance_block(group, rho, tau_study, tau_experiment):
    v = group.v.to_numpy(float)
    exp = group.Experiment_ID.astype(str).to_numpy()
    same_exp = (exp[:, None] == exp[None, :]).astype(float)
    sampling = np.diag(v) + rho * np.sqrt(np.outer(v, v)) * (same_exp - np.eye(len(v)))
    return sampling + tau_study * np.ones((len(v), len(v))) + tau_experiment * same_exp


def gls_components(data, rho, tau_study, tau_experiment):
    p = 2
    xtwx = np.zeros((p, p))
    xtwy = np.zeros(p)
    logdet = 0.0
    quad = 0.0
    pieces = []
    for study_id, group in data.groupby("Study_ID", sort=False):
        y = group.y.to_numpy(float)
        X = np.column_stack([np.ones(len(group)), group.warm.to_numpy(float)])
        cov = covariance_block(group, rho, tau_study, tau_experiment)
        try:
            cf = cho_factor(cov, lower=True, check_finite=False)
        except np.linalg.LinAlgError:
            return None
        wX = cho_solve(cf, X, check_finite=False)
        wy = cho_solve(cf, y, check_finite=False)
        xtwx += X.T @ wX
        xtwy += X.T @ wy
        logdet += 2 * np.log(np.diag(cf[0])).sum()
        quad += float(y @ wy)
        pieces.append((study_id, group, X, cf))
    try:
        beta = np.linalg.solve(xtwx, xtwy)
        sign, logdet_x = np.linalg.slogdet(xtwx)
    except np.linalg.LinAlgError:
        return None
    if sign <= 0:
        return None
    reml = logdet + logdet_x + quad - float(beta @ xtwy)
    return beta, xtwx, reml, pieces


def fit_meta_regression(data, rho=0.0, start=None):
    k = data.Study_ID.nunique()
    p = 2
    if k <= p:
        raise ValueError("At least three studies are required for this moderator model")
    scale = max(float(data.y.var()), float(data.v.median()), 0.01)

    def objective(x):
        ans = gls_components(data, rho, x[0], x[1])
        return 0.5 * ans[2] if ans is not None else 1e100

    starts = ([list(start)] if start is not None else
              [[0.0, 0.0], [scale / 4, scale / 4], [scale, 0.0], [0.0, scale]])
    fits = [minimize(objective, s, method="L-BFGS-B", bounds=[(0, None), (0, None)],
                     options={"maxiter": 500}) for s in starts]
    best = min(fits, key=lambda x: x.fun)
    if not np.isfinite(best.fun):
        raise RuntimeError("REML optimisation failed")
    tau_study, tau_exp = [max(float(x), 0.0) for x in best.x]
    beta, xtwx, _, pieces = gls_components(data, rho, tau_study, tau_exp)
    bread = np.linalg.inv(xtwx)
    meat = np.zeros((p, p))
    for _, group, X, cf in pieces:
        residual = group.y.to_numpy(float) - X @ beta
        score = X.T @ cho_solve(cf, residual, check_finite=False)
        meat += np.outer(score, score)
    robust = (k / (k - p)) * bread @ meat @ bread
    model = bread
    df = k - p
    critical = float(t.ppf(0.975, df))

    contrasts = {
        "Cold-origin": np.array([1.0, 0.0]),
        "Warm-origin": np.array([1.0, 1.0]),
        "Warm minus cold": np.array([0.0, 1.0]),
    }
    estimates = []
    for label, c in contrasts.items():
        estimate = float(c @ beta)
        se_model = float(np.sqrt(max(c @ model @ c, 0)))
        se_cluster = float(np.sqrt(max(c @ robust @ c, 0)))
        se = max(se_model, se_cluster)
        lo, hi = estimate - critical * se, estimate + critical * se
        p_value = float(2 * t.sf(abs(estimate / se), df)) if se > 0 else np.nan
        estimates.append(dict(contrast=label, estimate=estimate, SE_model=se_model,
                              SE_cluster=se_cluster, SE_used=se, CI_low=lo, CI_high=hi,
                              p_value=p_value, percent_change=100*np.expm1(estimate),
                              percent_low=100*np.expm1(lo), percent_high=100*np.expm1(hi)))
    return estimates, dict(studies=k, experiments=data.Experiment_ID.nunique(), effects=len(data),
                           rho=rho, tau2_study=tau_study, tau2_experiment=tau_exp,
                           df=df, convergence=bool(best.success))


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(INPUT)
    data = add_climatic_origin(data)
    data["y"] = pd.to_numeric(data.y, errors="raise")
    data["v"] = pd.to_numeric(data.v, errors="raise")

    feasibility = []
    for (driver, nutrient, origin), g in data.groupby(
            ["Driver_harmonised", "Nutrient_form", "Climatic_origin"], sort=True):
        feasibility.append(dict(Driver=driver, Nutrient=nutrient, Climatic_origin=origin,
                                studies=g.Study_ID.nunique(), experiments=g.Experiment_ID.nunique(),
                                effects=len(g)))
    pd.DataFrame(feasibility).to_csv(OUTPUT / "climatic_origin_feasibility.csv", index=False)

    result_rows = []
    loo_rows = []
    target_groups = [("Temperature", "Nitrate"), ("Irradiance", "Nitrate")]
    for driver, nutrient in target_groups:
        group = data.loc[data.Driver_harmonised.eq(driver) & data.Nutrient_form.eq(nutrient)].copy()
        for rho in (0.0, 0.5, 0.8):
            estimates, info = fit_meta_regression(group, rho=rho)
            cold = group.loc[group.Climatic_origin.eq("Cold-origin")]
            warm = group.loc[group.Climatic_origin.eq("Warm-origin")]
            for row in estimates:
                result_rows.append(dict(Driver=driver, Nutrient=nutrient,
                                        cold_studies=cold.Study_ID.nunique(),
                                        warm_studies=warm.Study_ID.nunique(), **info, **row))

        for omitted in sorted(group.Study_ID.unique()):
            reduced = group.loc[~group.Study_ID.eq(omitted)].copy()
            estimates, info = fit_meta_regression(reduced, rho=0.0)
            contrast = next(x for x in estimates if x["contrast"] == "Warm minus cold")
            loo_rows.append(dict(Driver=driver, Nutrient=nutrient, omitted_study=omitted,
                                 cold_studies=reduced.loc[reduced.Climatic_origin.eq("Cold-origin"), "Study_ID"].nunique(),
                                 warm_studies=reduced.loc[reduced.Climatic_origin.eq("Warm-origin"), "Study_ID"].nunique(),
                                 **info, **contrast))

    results = pd.DataFrame(result_rows)
    loo = pd.DataFrame(loo_rows)
    results.to_csv(OUTPUT / "climatic_origin_meta_regression_results.csv", index=False)
    loo.to_csv(OUTPUT / "climatic_origin_meta_regression_leave_one_out.csv", index=False)
    print(results.to_string(index=False))
    print("\nLeave-one-study-out warm-minus-cold ranges (rho=0):")
    print(loo.groupby(["Driver", "Nutrient"]).agg(
        estimate_min=("estimate", "min"), estimate_max=("estimate", "max"),
        ci_low_min=("CI_low", "min"), ci_low_max=("CI_low", "max"),
        ci_high_min=("CI_high", "min"), ci_high_max=("CI_high", "max"),
    ).to_string())


if __name__ == "__main__":
    main()
