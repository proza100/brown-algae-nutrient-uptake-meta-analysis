"""Run leave-one-study-out and irradiance-design sensitivity checks.

Execute from within the package after installing the listed requirements.
Outputs are written to the shared results directory; source data are read-only.
"""

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for _folder in (PROJECT_ROOT, PROJECT_ROOT / '02_meta_analysis'):
    if str(_folder) not in sys.path:
        sys.path.insert(0, str(_folder))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from meta_analysis_effects import prepare, fit_meta, OUTPUT


def leave_one_out(effects, rhos=(0.0,)):
    rows = []
    for driver in ("Temperature", "Irradiance"):
        parent = effects.loc[effects.Driver_harmonised.eq(driver)]
        for nutrient, group in [("All nutrients (exploratory)", parent), *parent.groupby("Nutrient_form")]:
            if group.Study_ID.nunique() < 3:
                continue
            for rho in rhos:
                full = fit_meta(group, rho)
                rows.append(dict(Driver=driver, Nutrient=nutrient,
                                 omitted_study="None (full model)", **full))
                for study_id in sorted(group.Study_ID.unique()):
                    reduced = group.loc[group.Study_ID.ne(study_id)]
                    if reduced.Study_ID.nunique() < 3:
                        continue
                    fit = fit_meta(reduced, rho,
                                   start=(full["tau2_study"], full["tau2_experiment"]))
                    label = group.loc[group.Study_ID.eq(study_id), "study_label"].iloc[0]
                    rows.append(dict(Driver=driver, Nutrient=nutrient,
                                     omitted_study=study_id, omitted_label=label,
                                     omitted_effects=int(group.Study_ID.eq(study_id).sum()), **fit))
    return pd.DataFrame(rows)


def light_design(effects, rhos=(0.0, 0.5, 0.8)):
    light = effects.loc[effects.Driver_harmonised.eq("Irradiance")]
    designs = [("Both direct designs", light),
               ("Intensity gradient only", light.loc[light.Environmental_contrast_type.eq("Irradiance intensity")]),
               ("Light versus darkness only", light.loc[light.Environmental_contrast_type.eq("Light versus darkness")])]
    rows = []
    for design, parent in designs:
        for nutrient, group in [("All nutrients (exploratory)", parent), *parent.groupby("Nutrient_form")]:
            k = group.Study_ID.nunique()
            if k < 3:
                rows.append(dict(Design=design, Nutrient=nutrient, studies=k,
                                 experiments=group.Experiment_ID.nunique(), effects=len(group),
                                 status="Not fitted: fewer than 3 studies"))
                continue
            for rho in rhos:
                fit = fit_meta(group, rho)
                rows.append(dict(Design=design, Nutrient=nutrient,
                                 status="Fitted", **fit))
    return pd.DataFrame(rows)


def plot_leave_one_out(table, driver, nutrient, rho=0.0):
    x = table.loc[table.Driver.eq(driver) & table.Nutrient.eq(nutrient) & table.rho.eq(rho)].copy()
    x = x.loc[x.omitted_study.ne("None (full model)")].sort_values("percent_change")
    if x.empty:
        return None
    fig, ax = plt.subplots(figsize=(10, max(4.5, 0.42 * len(x) + 1.8)))
    yy = np.arange(len(x))
    ax.errorbar(x.percent_change, yy,
                xerr=[x.percent_change - x.percent_low, x.percent_high - x.percent_change],
                fmt="o", markersize=6, capsize=3, color="#176B86")
    ax.axvline(0, color="#777777", linestyle="--", lw=1)
    full = table.loc[table.Driver.eq(driver) & table.Nutrient.eq(nutrient) &
                     table.rho.eq(rho) & table.omitted_study.eq("None (full model)")].iloc[0]
    ax.axvline(full.percent_change, color="#A34A16", lw=1.5, label="Full model")
    ax.set_yticks(yy, [f"Without {r.omitted_label} (m={int(r.omitted_effects)})" for r in x.itertuples()])
    ax.set_xlabel("Pooled change after omitting one study (%)", fontsize=12)
    ax.set_title(f"{driver}: {nutrient}; assumed ρ={rho:g}", fontsize=13)
    ax.tick_params(axis="both", labelsize=10)
    ax.legend(frameon=False, fontsize=10)
    fig.tight_layout()
    path = OUTPUT / f"leave_one_out_{driver.lower()}_{nutrient.lower().replace(' ', '_').replace('(', '').replace(')', '')}_rho{str(rho).replace('.', 'p')}.png"
    fig.savefig(path, dpi=600, bbox_inches="tight")
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    return path


if __name__ == "__main__":
    OUTPUT.mkdir(parents=True, exist_ok=True)
    effects, audit = prepare()
    loo = leave_one_out(effects)
    design = light_design(effects)
    loo.to_csv(OUTPUT / "sensitivity_leave_one_study_out.csv", index=False)
    design.to_csv(OUTPUT / "sensitivity_irradiance_design.csv", index=False)
    for driver in ("Temperature", "Irradiance"):
        plot_leave_one_out(loo, driver, "All nutrients (exploratory)")
        for nutrient in sorted(effects.loc[effects.Driver_harmonised.eq(driver), "Nutrient_form"].unique()):
            plot_leave_one_out(loo, driver, nutrient)
    print("LOO rows:", len(loo), "design rows:", len(design))
    print("Saved to", OUTPUT)
