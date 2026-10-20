# Table: public coordinate datasets in README.md.
from pathlib import Path
import argparse
import json

import prepare_data as data


def main():
    root = Path(__file__).resolve().parent
    expected = json.loads((root / "datasets.json").read_text())["datasets"]
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", choices=[row["name"] for row in expected])
    args = parser.parse_args()
    data.DATA = root / "rebuilt_data"
    data.RAW = data.DATA / "raw"
    data.RAW.mkdir(parents=True, exist_ok=True)
    checks = []
    for row in expected:
        if args.dataset and row["name"] != args.dataset:
            continue
        source = row["source"]
        data.obtain(source["url"], Path(source["file"]).name, source["sha256"])
        if row["name"] == "shorelines_europe":
            actual = data.shorelines()
        else:
            actual = data.roads(row["name"].split("_")[1].upper(), row["n"], row["center_lon_lat_degrees"])
        for field in ("n", "bytes", "sha256", "center_lon_lat_degrees", "radius_metres"):
            if actual[field] != row[field]:
                raise RuntimeError("Projected data mismatch: " + row["name"] + ": " + field)
        for field in ("sha256", "bytes"):
            if actual["source"][field] != source[field]:
                raise RuntimeError("Raw source mismatch: " + row["name"])
        checks.append({"dataset": row["name"], "n": row["n"], "source_sha256": source["sha256"],
                       "projected_sha256": row["sha256"], "status": "PASS"})
        print(json.dumps(checks[-1]), flush=True)
    (data.DATA / "verification.json").write_text(json.dumps(checks, indent=2) + "\n")


if __name__ == "__main__":
    main()
