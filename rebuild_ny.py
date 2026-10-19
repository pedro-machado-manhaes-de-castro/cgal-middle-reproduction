# Table 1: recover and project the DIMACS New York coordinate dataset.
from pathlib import Path
import json

import prepare_data


root = Path(__file__).resolve().parent
expected = next(x for x in json.loads((root / "datasets.json").read_text())["datasets"]
                if x["name"] == "roads_ny")
prepare_data.DATA = root / "rebuilt_ny"
prepare_data.RAW = prepare_data.DATA / "raw"
prepare_data.RAW.mkdir(parents=True, exist_ok=True)
prepare_data.obtain(expected["source"]["url"], "USA-road-d.NY.co.gz", expected["source"]["sha256"])
actual = prepare_data.roads("NY", 264346, [-74.0, 40.8])
if actual["source"]["sha256"] != expected["source"]["sha256"] or actual["n"] != expected["n"]:
    raise RuntimeError("Source or record count differs from the recorded input")
exact = actual["sha256"] == expected["sha256"]
result = {"raw_source_sha256_verified": True, "n": actual["n"],
          "projected_sha256": actual["sha256"], "projected_bytes_match_recorded": exact,
          "input_for_benchmark": "data/roads_ny.xy"}
(prepare_data.DATA / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
if not exact:
    raise RuntimeError("Projection rounding differs on this platform. Use the bundled, hash-verified data/roads_ny.xy to reproduce the exact benchmark input.")
