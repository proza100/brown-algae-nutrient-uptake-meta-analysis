# Packaging checks â€” 4 November 2026

- All seven Python files parsed successfully.
- The command-line primary workflow completed using Python 3.13 and the dependency versions in requirements-tested.txt. Openpyxl was loaded from the bundled local library because the system environment did not include it; a clean pip installation was not tested (network unavailable).
- Exactly 358 selected Effect_IDs matched the existing analytical CSV (298 temperature, 60 irradiance). Oriented lnRR values and sampling variances matched exactly.
- All 24 primary/rho-sensitivity models were regenerated. Maximum absolute differences from saved results: lnRR 2.75e-6; confidence limits 1.48e-5; percentage change 0.000330 percentage points. Convergence flags were unchanged. These small numerical differences are consistent with numerical environment/optimization variation; bitwise identity is not claimed.
- Heterogeneity export completed.
- The environmental R workflow completed with R 4.5.2, producing 39 fitted-model rows, 15 feasibility rows and 78 leave-one-study-out rows. R printed locale startup warnings but completed successfully.
- Computational function bodies in forest plotting, Python sensitivity scripts and climatic-origin regression match the current working-project scripts. The core prepare function differs only in its portable worksheet reference; the resulting selected data were checked as above.
- The Python leave-one-study-out, one-effect selection, climatic-origin and full plotting stages were not rerun end-to-end during this packaging check. They are retained from the existing portable package, with unchanged computational functions. This is a packaging validation, not a new scientific audit.
- Input files, original scripts and original results were not modified. Validation outputs were written to a separate temporary directory and are not distributed in this ZIP.

Version 2 includes the two environmental context CSV inputs in data/; SD1 is downloaded separately. No previously computed results are included. INPUTS.sha256 identifies the exact inputs used in this check; context_columns.json lists their headers without including data.

## Version 2 packaging update

The two included CSVs were checked byte-for-byte against the validated original inputs. Default file paths now use data/Supplementary_Data_SD1.xlsx and the included CSVs. Download instructions and missing-file guidance were added. No statistical functions changed; models were not rerun for this file-location update.
