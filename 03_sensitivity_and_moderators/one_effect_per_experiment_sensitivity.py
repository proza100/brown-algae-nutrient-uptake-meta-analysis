"""Deterministic one-effect-per-experiment-and-nutrient sensitivity analysis.

Selection never uses lnRR, statistical significance, or its sampling variance.
See selection_rule() for the fully specified hierarchy.
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

from meta_analysis_effects import OUTPUT, prepare, fit_meta


def selection_rule(effects, species_tiebreak="first"):
    """Select one effect per Driver × Experiment_ID × Nutrient_form.

    Priority: original-study primary response; positive exposure change when
    available; median *unique* absolute exposure contrast within the group;
    latest measured time point; alphabetically first/last species; lowest
    numeric Effect_ID. Missing exposure differences rank after known ones.
    This is a sensitivity rule, not a claim that the chosen contrast was the
    authors' primary planned comparison.
    """
    if species_tiebreak not in ("first", "last"):
        raise ValueError("species_tiebreak must be first or last")
    selected = []
    keys = ["Driver_harmonised", "Experiment_ID", "Nutrient_form"]
    for _, group in effects.groupby(keys, sort=True):
        g = group.copy()
        role = g.Outcome_role_in_original_study.fillna("").str.strip().str.lower()
        g["_role_rank"] = np.where(role.eq("primary_response"), 0, 1)
        direction = g.Environmental_change_direction.fillna("").str.lower()
        g["_direction_rank"] = np.where(direction.str.startswith("higher"), 0, 1)
        difference = pd.to_numeric(g.Exposure_difference, errors="coerce").abs()
        unique_magnitudes = np.sort(difference.dropna().unique())
        median_magnitude = np.median(unique_magnitudes) if len(unique_magnitudes) else np.nan
        g["_contrast_rank"] = (difference - median_magnitude).abs().fillna(np.inf)
        g["_latest_rank"] = -pd.to_numeric(g.Timepoint_value, errors="coerce").fillna(-np.inf)
        g["_species_rank"] = g.Species.fillna("").str.casefold()
        g["_effect_id_rank"] = pd.to_numeric(g.Effect_ID.astype(str).str.extract(r"(\d+)")[0],
                                              errors="coerce").fillna(np.inf)
        g = g.sort_values(["_role_rank", "_direction_rank", "_contrast_rank",
                           "_latest_rank", "_species_rank", "_effect_id_rank"],
                          ascending=[True, True, True, True,
                                     species_tiebreak == "first", True], kind="stable")
        selected.append(g.iloc[0])
    out = pd.DataFrame(selected).drop(columns=[c for c in g.columns if c.startswith("_")])
    if out.duplicated(keys).any():
        raise AssertionError("Selection did not produce unique experiment × nutrient rows")
    return out.sort_values(keys).reset_index(drop=True)


def compare_models(full, selected, rho_full=0.0):
    rows = []
    for driver in ("Temperature", "Irradiance"):
        a = full.loc[full.Driver_harmonised.eq(driver)]
        b = selected.loc[selected.Driver_harmonised.eq(driver)]
        for nutrient in ["All nutrients (exploratory)", *sorted(a.Nutrient_form.unique())]:
            af = a if nutrient == "All nutrients (exploratory)" else a.loc[a.Nutrient_form.eq(nutrient)]
            bs = b if nutrient == "All nutrients (exploratory)" else b.loc[b.Nutrient_form.eq(nutrient)]
            if af.Study_ID.nunique() < 3 or bs.Study_ID.nunique() < 3:
                continue
            for name, data, rho in [("All selected effects", af, rho_full),
                                    ("One per experiment × nutrient", bs, 0.0)]:
                fit = fit_meta(data, rho)
                rows.append(dict(Driver=driver, Nutrient=nutrient, Analysis=name, **fit))
    return pd.DataFrame(rows)


def comparison_plot(results, driver, outdir=OUTPUT):
    df = results.loc[results.Driver.eq(driver)].copy()
    df["order"] = df.Nutrient.map({"All nutrients (exploratory)": 0,
                                    "Nitrate": 1, "Ammonium": 2, "Phosphate": 3})
    df = df.sort_values(["order", "Analysis"])
    nutrients = df.Nutrient.drop_duplicates().tolist()
    fig, ax = plt.subplots(figsize=(9, max(4.2, 0.85 * len(nutrients) + 1.4)))
    for i, nutrient in enumerate(nutrients):
        subset = df.loc[df.Nutrient.eq(nutrient)]
        for row in subset.itertuples():
            short = row.Analysis.startswith("One per")
            y = i + (0.13 if short else -0.13)
            color = "#B0522E" if short else "#1A6B83"
            ax.errorbar(row.percent_change, y,
                        xerr=[[row.percent_change - row.percent_low],
                              [row.percent_high - row.percent_change]],
                        fmt="s" if short else "o", color=color,
                        capsize=3, ms=5, lw=1.2)
    ax.axvline(0, color="#666666", ls="--", lw=1)
    ax.set_yticks(range(len(nutrients)), nutrients)
    ax.invert_yaxis()
    ax.set_xlabel("Estimated change in nutrient uptake (%) with higher exposure")
    ax.set_title(f"{driver}: all effects versus one per experiment and nutrient")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0], [0], marker="o", color="#1A6B83", lw=0,
                              label="All selected effects"),
                       Line2D([0], [0], marker="s", color="#B0522E", lw=0,
                              label="One per experiment × nutrient")],
              frameon=False, loc="best")
    ax.grid(axis="x", color="#E4E4E4", linewidth=0.6)
    fig.tight_layout()
    path = Path(outdir) / f"sensitivity_one_effect_per_experiment_{driver.lower()}.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


if __name__ == "__main__":
    OUTPUT.mkdir(parents=True, exist_ok=True)
    effects, _ = prepare()
    chosen = selection_rule(effects)
    alternate = selection_rule(effects, species_tiebreak="last")
    chosen.to_csv(OUTPUT / "one_effect_per_experiment_selection.csv", index=False)
    alternate_ids = alternate[["Driver_harmonised", "Experiment_ID", "Nutrient_form", "Effect_ID"]]
    chosen_ids = chosen[["Driver_harmonised", "Experiment_ID", "Nutrient_form", "Effect_ID"]]
    differences = chosen_ids.merge(alternate_ids,
                                   on=["Driver_harmonised", "Experiment_ID", "Nutrient_form"],
                                   suffixes=("_first_species", "_last_species"))
    differences = differences.loc[differences.Effect_ID_first_species.ne(
        differences.Effect_ID_last_species)]
    differences.to_csv(OUTPUT / "one_effect_selection_species_tiebreak_audit.csv", index=False)
    results = compare_models(effects, chosen)
    results.to_csv(OUTPUT / "sensitivity_one_effect_per_experiment_results.csv", index=False)
    alternate_results = compare_models(effects, alternate)
    alternate_results.to_csv(OUTPUT / "sensitivity_one_effect_per_experiment_alternate_species.csv",
                             index=False)
    for driver in ("Temperature", "Irradiance"):
        print(comparison_plot(results, driver))
    print("Selected", len(chosen), "of", len(effects), "effects")
    print("Species tie-break differences:", len(differences))
