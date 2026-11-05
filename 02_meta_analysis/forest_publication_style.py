"""Publication-style forest plots with experiment rows and nutrient subgroups.

Layout is inspired by the supplied reference figure; statistics remain the
project's lnRR models, not the reference figure's NMD analysis.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Polygon
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import chi2

from meta_analysis_effects import OUTPUT, prepare, _display_points, _cov_block


COLOURS = {"Nitrate": "#237A3B", "Ammonium": "#B44D2A", "Phosphate": "#263F99"}
NUTRIENTS = ("Nitrate", "Ammonium", "Phosphate")


def _fmt(value):
    return f"{value:+.2f}" if np.isfinite(value) else "—"


def _diamond(ax, y, mean, low, high, color, height=0.25):
    if not np.all(np.isfinite([mean, low, high])):
        return
    ax.add_patch(Polygon([(low, y), (mean, y + height), (high, y),
                          (mean, y - height)], closed=True,
                         facecolor=color, edgecolor=color, linewidth=0.8, zorder=5))


def heterogeneity_diagnostics(data, rho, tau2_study, tau2_experiment):
    """Generalized fixed-effect Q and approximate multilevel I².

    The assumed sampling covariance V contains rho within experiments. Typical
    sampling variance is (n-1)/tr(P), where P is the GLS residual projector.
    Q/p are diagnostics conditional on V, not independent replication tests.
    """
    total_w1 = total_wy = sum_ywy = trace_w = sum_w1sq = 0.0
    for _, group in data.groupby("Study_ID", sort=False):
        v = _cov_block(group, rho, 0., 0.)
        cf = cho_factor(v, lower=True, check_finite=False)
        one = np.ones(len(group))
        y = group.y.to_numpy(float)
        w1 = cho_solve(cf, one, check_finite=False)
        wy = cho_solve(cf, y, check_finite=False)
        total_w1 += float(one @ w1)
        total_wy += float(one @ wy)
        sum_ywy += float(y @ wy)
        trace_w += float(np.trace(cho_solve(cf, np.eye(len(group)), check_finite=False)))
        sum_w1sq += float(w1 @ w1)
    n = len(data)
    q = max(0.0, sum_ywy - total_wy ** 2 / total_w1)
    df = n - 1
    trace_p = trace_w - sum_w1sq / total_w1
    typical_v = df / trace_p if trace_p > 0 else np.nan
    tau2_total = tau2_study + tau2_experiment
    i2 = 100 * tau2_total / (tau2_total + typical_v) if typical_v > 0 else np.nan
    return dict(Q=q, df=df, p=chi2.sf(q, df), typical_v=typical_v,
                tau2_total=tau2_total, I2_approx=i2)


def _p_text(p):
    return "<0.001" if p < 0.001 else f"={p:.3f}"


def build_figure(effects, results, driver, rho=0.0, outdir=OUTPUT):
    """One row per experiment × nutrient; model diamonds from effect-level REML."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    subset = effects.loc[effects.Driver_harmonised.eq(driver)]
    rows = []
    for nutrient in NUTRIENTS:
        nutrient_effects = subset.loc[subset.Nutrient_form.eq(nutrient)]
        if nutrient_effects.empty:
            continue
        rows.append({"kind": "header", "nutrient": nutrient})
        points = _display_points(nutrient_effects, rho, "Experiment")
        points = points.sort_values(["lnRR", "label", "group_id"], ascending=True)
        for p in points.itertuples(index=False):
            ex = nutrient_effects.loc[nutrient_effects.Experiment_ID.eq(p.group_id)]
            species = sorted(ex.Species.dropna().astype(str).unique())
            species_text = "\n".join(species)
            rows.append(dict(kind="experiment", nutrient=nutrient,
                             study=ex.study_label.iloc[0], experiment=p.group_id.split("_")[-1],
                             species=species_text, m=p.effects, mean=p.lnRR,
                             low=p.lnRR - 1.96 * p.SE, high=p.lnRR + 1.96 * p.SE))
        fit = results.loc[results.Driver.eq(driver) & results.Nutrient.eq(nutrient) &
                          np.isclose(results.rho, rho)]
        if len(fit):
            r = fit.iloc[0]
            rows.append(dict(kind="subgroup", nutrient=nutrient, mean=r.lnRR,
                             low=r.CI_low, high=r.CI_high,
                             studies=int(r.studies), experiments=int(r.experiments)))
            rows.append(dict(kind="heterogeneity", nutrient=nutrient,
                             tau2_study=r.tau2_study,
                             tau2_experiment=r.tau2_experiment,
                             **heterogeneity_diagnostics(nutrient_effects, rho,
                                                         r.tau2_study, r.tau2_experiment)))
            rows.append(dict(kind="q_diagnostic", nutrient=nutrient,
                             **heterogeneity_diagnostics(nutrient_effects, rho,
                                                         r.tau2_study, r.tau2_experiment)))
        else:
            rows.append(dict(kind="no_model", nutrient=nutrient,
                             studies=nutrient_effects.Study_ID.nunique()))
        rows.append({"kind": "space", "nutrient": nutrient})
    overall = results.loc[results.Driver.eq(driver) &
                          results.Nutrient.eq("All nutrients (exploratory)") &
                          np.isclose(results.rho, rho)]
    if len(overall):
        r = overall.iloc[0]
        rows.append(dict(kind="overall", nutrient="All nutrients", mean=r.lnRR,
                         low=r.CI_low, high=r.CI_high, studies=int(r.studies)))
        rows.append(dict(kind="heterogeneity", nutrient="All nutrients",
                         tau2_study=r.tau2_study,
                         tau2_experiment=r.tau2_experiment,
                         **heterogeneity_diagnostics(subset, rho,
                                                     r.tau2_study, r.tau2_experiment)))
        rows.append(dict(kind="q_diagnostic", nutrient="All nutrients",
                         **heterogeneity_diagnostics(subset, rho,
                                                     r.tau2_study, r.tau2_experiment)))

    values = [v for row in rows for v in (row.get("low"), row.get("high"))
              if v is not None and np.isfinite(v)]
    left, right = min(values + [0]), max(values + [0])
    span = max(right - left, 0.6)
    left -= span * 0.06
    right += span * 0.06
    n = len(rows)
    fig = plt.figure(figsize=(19.5, max(8.0, 0.43 * n + 2.1)), dpi=160,
                     facecolor="white")
    grid = GridSpec(1, 2, figure=fig, width_ratios=[1.6, 1.1], wspace=0.015,
                    left=0.035, right=0.985, top=0.91, bottom=0.11)
    tab = fig.add_subplot(grid[0, 0])
    ax = fig.add_subplot(grid[0, 1], sharey=tab)
    tab.set_xlim(0, 1)
    tab.set_ylim(n + 0.6, -0.8)
    tab.axis("off")
    ax.set_xlim(left, right)
    ax.set_ylim(n + 0.6, -0.8)
    ax.set_yticks([])
    for spine in ("left", "right", "top"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#555555")
    ax.axvline(0, color="#888888", linestyle=":", linewidth=1.2, zorder=0)
    ax.set_xlabel("Log response ratio (lnRR): higher exposure versus lower exposure", fontsize=16)
    ax.tick_params(axis="x", labelsize=10)
    columns = [(0.00, "Study"), (0.29, "Experiment"), (0.43, "Species"),
               (0.72, "Effects (m)"), (0.82, "lnRR [95% CI]")]
    for xpos, label in columns:
        tab.text(xpos, -0.9, label, fontsize=14, fontweight="bold", va="bottom")
    ax.text((left + right) / 2, -0.9, "Effect and 95% CI", fontsize=14,
            fontweight="bold", ha="center", va="bottom")
    tab.plot([0, 1], [-0.45, -0.45], color="#555555", lw=0.7, clip_on=False)
    ax.plot([left, right], [-0.45, -0.45], color="#555555", lw=0.7, clip_on=False)
    for y, row in enumerate(rows):
        kind = row["kind"]
        nut = row["nutrient"]
        col = COLOURS.get(nut, "#202020")
        if kind == "header":
            tab.text(0, y, nut, fontsize=14, fontweight="bold", va="center")
            tab.plot([0, 1], [y + 0.42, y + 0.42], color="#DDDDDD", lw=0.5)
            ax.plot([left, right], [y + 0.42, y + 0.42], color="#DDDDDD", lw=0.5)
        elif kind == "experiment":
            tab.text(0, y, row["study"], fontsize=14, va="center")
            tab.text(0.29, y, row["experiment"], fontsize=14, va="center")
            tab.text(0.43, y, row["species"], fontsize=14, va="center",
                     style="italic", linespacing=0.95)
            tab.text(0.765, y, str(row["m"]), fontsize=14, ha="center", va="center")
            tab.text(0.82, y, f"{_fmt(row['mean'])} [{_fmt(row['low'])}, {_fmt(row['high'])}]",
                     fontsize=14, va="center")
            ax.plot([row["low"], row["high"]], [y, y], color="#505050", lw=0.9, zorder=2)
            ax.plot(row["mean"], y, marker="s", color=col, markersize=6.2, zorder=3)
        elif kind == "subgroup":
            tab.text(0, y, f"{nut} pooled ({row['studies']} studies; {row['experiments']} experiments)",
                     fontsize=14, fontweight="bold", va="center")
            tab.text(0.82, y, f"{_fmt(row['mean'])} [{_fmt(row['low'])}, {_fmt(row['high'])}]",
                     fontsize=14, fontweight="bold", va="center")
            _diamond(ax, y, row["mean"], row["low"], row["high"], col)
        elif kind == "no_model":
            tab.text(0, y, f"{nut}: no pooled estimate ({row['studies']} studies)",
                     fontsize=14, color="#555555", va="center")
        elif kind == "overall":
            tab.text(0, y, f"All nutrients pooled ({row['studies']} studies; exploratory)",
                     fontsize=14, fontweight="bold", va="center")
            tab.text(0.82, y, f"{_fmt(row['mean'])} [{_fmt(row['low'])}, {_fmt(row['high'])}]",
                     fontsize=14, fontweight="bold", va="center")
            _diamond(ax, y, row["mean"], row["low"], row["high"], "#222222")
        elif kind == "heterogeneity":
            summary = (f"Heterogeneity: τ² study={row['tau2_study']:.3f}; "
                       f"τ² experiment={row['tau2_experiment']:.3f}; "
                       f"I² total≈{row['I2_approx']:.1f}%")
            tab.text(0, y, summary, fontsize=12, color="#555555", va="center")
        elif kind == "q_diagnostic":
            diagnostic = f"Q({row['df']})={row['Q']:.1f}, p{_p_text(row['p'])} (diagnostic; assumes displayed ρ)"
            tab.text(0, y, diagnostic, fontsize=12, color="#555555", va="center")
    fig.suptitle(f"{driver} and nutrient uptake", fontsize=16, fontweight="bold", y=0.972)
    fig.text(0.035, 0.045,
             f"Experiment points are descriptive GLS summaries (m = number of effects); diamonds use the effect-level multilevel REML models. Assumed sampling ρ = {rho:g}.",
             fontsize=12, color="#404040")
    fig.text(0.035, 0.026,
             "I² is approximate; Q/p depend on assumed sampling covariance. CI for experiment points assumes the displayed ρ. Overall cross-nutrient estimate is exploratory.",
             fontsize=12, color="#404040")
    stem = f"forest_{driver.lower()}_publication_style_rho{str(rho).replace('.', 'p')}"
    png = outdir / f"{stem}.png"
    tiff = outdir / f"{stem}_300dpi.tiff"
    fig.savefig(png, dpi=220, bbox_inches="tight", facecolor="white")
    if rho == 0:
        fig.savefig(tiff, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return png, (tiff if rho == 0 else None)


def export_heterogeneity_table(effects, results, outdir=OUTPUT):
    records = []
    for r in results.itertuples(index=False):
        data = effects.loc[effects.Driver_harmonised.eq(r.Driver)]
        if r.Nutrient != "All nutrients (exploratory)":
            data = data.loc[data.Nutrient_form.eq(r.Nutrient)]
        diagnostics = heterogeneity_diagnostics(data, r.rho,
                                                r.tau2_study, r.tau2_experiment)
        records.append(dict(Driver=r.Driver, Nutrient=r.Nutrient, rho=r.rho,
                            studies=r.studies, experiments=r.experiments,
                            effects=r.effects, tau2_study=r.tau2_study,
                            tau2_experiment=r.tau2_experiment, **diagnostics))
    output = Path(outdir) / "meta_analysis_heterogeneity.csv"
    pd.DataFrame(records).to_csv(output, index=False)
    return output


if __name__ == "__main__":
    effects, _ = prepare()
    results = pd.read_csv(OUTPUT / "meta_analysis_results.csv", sep=None, engine="python")
    print(export_heterogeneity_table(effects, results))
    for driver in ("Temperature", "Irradiance"):
        for rho in (0.0, 0.5, 0.8):
            print(build_figure(effects, results, driver, rho))
