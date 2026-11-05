# Brown macroalgal nutrient uptake: meta-analysis code

Analysis code and environmental-input publication package, 4 November 2026. Includes primary multilevel models, heterogeneity diagnostics, rho sensitivity (0, 0.5, 0.8), leave-one-study-out analyses, irradiance-design sensitivity, one effect per experiment × nutrient within each driver, climatic-origin meta-regressions, environmental meta-regressions, and forest plots. Descriptive-review figures, extraction tools, Excel builders, exploratory measurement-method subgroups, notebooks, the separately published SD1 workbook and precomputed results are excluded.

## Input data and obtaining SD1

The `data/` directory includes `temperature_effects_detail.csv` and `irradiance_effects_detail.csv`. These contain the numerical exposure/acclimation metadata used in environmental regressions and are joined by Effect_ID, Study_ID and Experiment_ID. They are read automatically; no separate download is required.

The analytical Excel workbook is published separately as **Supplementary Data SD1** and is not duplicated in this ZIP:

1. Open the article's publication page and its Supporting Information / Supplementary Material section (or the associated data repository, if that is where the journal provides SD1).
2. Download **Supplementary Data SD1**, filename `Supplementary_Data_SD1.xlsx`. Download the actual Excel attachment rather than the supplementary Methods document or SD4 results.
3. Place it in this package's `data/` directory, alongside the two CSV files. If the publisher changes the downloaded filename, rename it to `Supplementary_Data_SD1.xlsx` or supply its actual path with `--data`.
4. Keep worksheet and column names unchanged. The default input sheet is `All_outcomes`; another name can be supplied with `--sheet`.

No publication URL or DOI was available when this package was prepared; no download address is assumed. A direct SD1 link can be added to this README once the publication or repository landing page is available.

The workbook contains 566 outcome records. Original selection/orientation rules recover 358 effects (298 temperature, 60 irradiance). This workflow begins with extracted effects and does not re-extract data from papers. INPUTS.sha256 identifies the inputs used for the packaging validation; the SD1 hash refers to the previously named workbook with unchanged contents. A hash mismatch after substantive data revisions requires checking the input version.

## Installation and execution

Python 3.13 was used for the packaging check. Install `requirements.txt`, or use `requirements-tested.txt` for the exact installed dependency versions. The latter describes the packaging environment, not a claim about the original analysis environment. Environmental regressions require R (base R only) and Rscript on PATH.

```sh
python -m pip install -r requirements-tested.txt
# After downloading SD1 into data/:
python run_analysis.py --stage primary
python run_analysis.py --stage sensitivity
python run_analysis.py --stage climate
python run_analysis.py --stage forests
python run_analysis.py --stage all
# Alternatively, use SD1 from another folder:
python run_analysis.py --data "/path/to/Supplementary_Data_SD1.xlsx" --stage all
```

Quote paths containing spaces. Commands work from any working directory when the script path is specified. Each stage regenerates primary results first; files are written beneath `results/`, or the directory supplied with `--output`. Inputs are read-only. Use a fresh output directory to retain previous runs. No proprietary software, Jupyter or OpenAI runtime is required.

## Scientific conventions and limitations

Positive lnRR indicates greater uptake at higher exposure after orientation. Covariance assumptions, fitting algorithms, starting values and inference remain as in the retained analysis scripts. Estimates are not standardised responses to a fixed temperature or irradiance increase. Overall cross-nutrient estimates are exploratory; the irradiance overall model retains nitrite and urea. Separate study and experiment variance components are not identifiable when each study contributes only one experiment; interpret their sum in those groups. Optimizer convergence flags must be inspected; packaging does not resolve previously recorded convergence issues. Environmental moderator slopes are per 1 °C or per 100 µmol photons m−2 s−1, and their leave-one-study-out analyses use rho=0.

## Provenance

The five Python analysis modules are copied from the portable September 2026 code package. The environmental regression script is added from the existing environmental-contrast analysis; only its path configuration is changed. `project_config.py` and `run_analysis.py` provide portable execution. Original project data, scripts and results are not modified. `MANIFEST.sha256` records package file hashes. See `VALIDATION.md` for actual checks and their limits.
