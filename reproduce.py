# Tables 1-3 and Figures 1-3: reproduce the Delaunay policy comparison.
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys

root = Path(__file__).resolve().parent
os.chdir(root)
base = Path(os.environ["CGAL_BASE"]).resolve()
manifest = json.loads((root / "measurement_experiment.json").read_text())
original = (root / "headers/original/CGAL/Hilbert_sort_middle_2.h").read_bytes()
base_header = (base / "Spatial_sorting/include/CGAL/Hilbert_sort_middle_2.h").read_bytes()
if base_header.replace(b"\r\n", b"\n") != original.replace(b"\r\n", b"\n"):
    raise RuntimeError("CGAL middle header differs from the pinned original source")
head = subprocess.check_output(["git", "-C", str(base), "rev-parse", "HEAD"], text=True).strip()
if head != manifest["base_commit"]:
    raise RuntimeError("CGAL checkout is not at the pinned commit")
if subprocess.run(["git", "-C", str(base), "diff", "--quiet", "HEAD", "--"]).returncode:
    raise RuntimeError("CGAL checkout contains tracked modifications")
manifest["base_header_sha256"] = hashlib.sha256(base_header).hexdigest()
manifest["cgal_base"] = str(base)
manifest["line_ending_adaptation"] = "The original source was compared byte-for-byte after CRLF-to-LF normalization. Only the expected local checkout hash was adapted."
(root / "experiment.json").write_text(json.dumps(manifest, indent=2) + "\n")
subprocess.run([sys.executable, "run_campaign.py", "--phase", "prepare"], check=True)
expected = json.loads((root / "datasets.json").read_text())["datasets"]
actual = {x["name"]: x for x in json.loads((root / "data/datasets.json").read_text())["datasets"]}
for dataset in expected:
    converted = actual[dataset["name"]]
    for field in ["n", "sha256"]:
        if converted[field] != dataset[field]:
            raise RuntimeError("Converted public data differs from the recorded inputs")
    if converted["source"]["sha256"] != dataset["source"]["sha256"]:
        raise RuntimeError("Downloaded public source differs from the recorded release")
subprocess.run([sys.executable, "run_campaign.py", "--phase", "run"], check=True)
subprocess.run([sys.executable, "analyze.py"], check=True)
