# Data required for execution

Included: `temperature_effects_detail.csv` and `irradiance_effects_detail.csv` (environmental moderator inputs).

Download `Supplementary_Data_SD1.xlsx` from the article's Supporting Information / Supplementary Material entry **Supplementary Data SD1**, or the associated repository if provided by the journal. Place the downloaded workbook in this directory. No direct link was available at packaging; consult the article's publication page. Do not substitute SD4, which contains model results rather than the input workbook.

Expected directory contents after downloading:

```
data/
  README.md
  Supplementary_Data_SD1.xlsx
  temperature_effects_detail.csv
  irradiance_effects_detail.csv
```

From the package root, run `python run_analysis.py --stage all`. Python and R must be installed as described in the main README. The CSV paths are resolved relative to the package automatically. SD1 may instead be supplied from elsewhere using `--data`. No input file is modified.
