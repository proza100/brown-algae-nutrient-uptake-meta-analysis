"""Portable input/output configuration; input data are distributed separately."""
from pathlib import Path
import os
PROJECT_ROOT = Path(__file__).resolve().parent
DATA_FILE = Path(os.environ.get("ALGAE_DATA", PROJECT_ROOT / "data" / "Supplementary_Data_SD1.xlsx")).resolve()
DATA_SHEET = os.environ.get("ALGAE_SHEET", "All_outcomes")
RESULTS_DIR = Path(os.environ.get("ALGAE_RESULTS", PROJECT_ROOT / "results")).resolve()
