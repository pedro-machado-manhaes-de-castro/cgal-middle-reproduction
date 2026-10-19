# Tables 2-3 and Figures 1-3: regenerate summaries and figures from recorded Delaunay trials.
from pathlib import Path
import subprocess
import sys

import analyze


root = Path(__file__).resolve().parent
analyze.RESULTS = root / "recorded_results"
analyze.OUTPUT = root / "analysis"
analyze.OUTPUT.mkdir(exist_ok=True)
analyze.main()
subprocess.run([sys.executable, str(root / "plot_figures.py"),
                "--analysis", str(analyze.OUTPUT), "--output", str(root / "figures")], check=True)
