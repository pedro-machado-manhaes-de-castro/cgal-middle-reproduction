# Tables 2-3: verify the summaries of all 2730 recorded Delaunay constructions.
from pathlib import Path
import csv
import hashlib
import json
import math

import analyze


root = Path(__file__).resolve().parent
analyze.RESULTS = root / "recorded_results"
analyze.OUTPUT = root / "rechecked_analysis"
analyze.OUTPUT.mkdir(exist_ok=True)
analyze.main()
files = ("seed_measurements.csv", "timings_and_paths.csv", "paired_ratios.csv", "largest_case_ratios.csv")
checks = {}
for name in files:
    original = root / "analysis" / name
    rebuilt = analyze.OUTPUT / name
    expected = list(csv.DictReader(original.open(newline="", encoding="utf-8")))
    actual = list(csv.DictReader(rebuilt.open(newline="", encoding="utf-8")))
    if len(actual) != len(expected):
        raise RuntimeError("Summary row count differs: " + name)
    for a, b in zip(actual, expected):
        if a.keys() != b.keys():
            raise RuntimeError("Summary columns differ: " + name)
        for key in a:
            if a[key] != b[key]:
                try:
                    same = math.isclose(float(a[key]), float(b[key]), rel_tol=1e-13, abs_tol=1e-15)
                except ValueError:
                    same = False
                if not same:
                    raise RuntimeError("Summary value differs: " + name + ": " + key)
    checks[name] = {"rows": len(actual), "numerically_equal": True,
                   "byte_identical": original.read_bytes() == rebuilt.read_bytes(),
                   "recorded_sha256": hashlib.sha256(original.read_bytes()).hexdigest()}
result = {"status": "PASS", "checks": checks, "measurement_reruns": 0}
(analyze.OUTPUT / "recorded_verification.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
