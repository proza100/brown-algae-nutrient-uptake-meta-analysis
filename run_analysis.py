"""Run the published models without Jupyter. See README.md for external inputs."""
import argparse, os, sys, runpy, subprocess
from pathlib import Path
root=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("--data", type=Path, default=root/"data"/"Supplementary_Data_SD1.xlsx")
p.add_argument("--sheet", default="All_outcomes")
p.add_argument("--output", type=Path, default=Path(__file__).resolve().parent/"results")
p.add_argument("--stage", choices=["primary","sensitivity","climate","forests","environment","all"],default="primary")
p.add_argument("--temperature-context",type=Path,default=root/"data"/"temperature_effects_detail.csv")
p.add_argument("--irradiance-context",type=Path,default=root/"data"/"irradiance_effects_detail.csv")
p.add_argument("--rscript",default="Rscript")
a=p.parse_args()
if not a.data.is_file(): p.error("Download Supplementary_Data_SD1.xlsx from the article Supporting Information (Supplementary Data SD1), place it in data/, or supply --data. See data/README.md.")
if a.stage in ("environment","all") and (not a.temperature_context.is_file() or not a.irradiance_context.is_file()): p.error("Environmental regressions require both context CSV files")
os.environ["ALGAE_DATA"]=str(a.data.resolve());os.environ["ALGAE_SHEET"]=a.sheet
os.environ["ALGAE_RESULTS"]=str(a.output.resolve());os.environ.setdefault("MPLBACKEND","Agg")
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/"02_meta_analysis"),str(root/"03_sensitivity_and_moderators")]
from meta_analysis_effects import prepare,run_models,OUTPUT
from forest_publication_style import export_heterogeneity_table
effects,audit=prepare();OUTPUT.mkdir(parents=True,exist_ok=True)
if len(effects)!=358 or effects.Driver_harmonised.value_counts().to_dict()!={"Temperature":298,"Irradiance":60}: raise ValueError("Input does not reproduce the published 358-effect selection")
effects.to_csv(OUTPUT/"analysis_input_effects.csv",index=False)
audit.to_csv(OUTPUT/"effect_selection_audit.csv",index=False)
results=run_models(effects,rhos=(0.0,0.5,0.8))
results.to_csv(OUTPUT/"meta_analysis_results.csv",index=False)
export_heterogeneity_table(effects,results)
if a.stage in ("sensitivity","all"):
 for f in ["run_meta_sensitivity.py","one_effect_per_experiment_sensitivity.py"]: runpy.run_path(str(root/"03_sensitivity_and_moderators"/f),run_name="__main__")
if a.stage in ("climate","all"): runpy.run_path(str(root/"03_sensitivity_and_moderators/climatic_origin_meta_regression.py"),run_name="__main__")
if a.stage in ("forests","all"): runpy.run_path(str(root/"02_meta_analysis/forest_publication_style.py"),run_name="__main__")
if a.stage in ("environment","all"):
 subprocess.run([a.rscript,str(root/"03_sensitivity_and_moderators/environmental_meta_regressions.R"),str(OUTPUT/"analysis_input_effects.csv"),str(a.temperature_context.resolve()),str(a.irradiance_context.resolve()),str(a.output.resolve()/"environmental_meta_regressions")],check=True)
print("Results:",a.output.resolve())
