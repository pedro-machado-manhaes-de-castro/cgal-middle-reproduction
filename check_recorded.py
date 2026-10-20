# Table: all 46 corrected/original total-time comparisons in README.md.
from collections import Counter
from pathlib import Path
import csv
import json
import math

import analyze


ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def samples(directory):
    return [json.loads(line) for line in (ROOT / directory / "samples.jsonl").read_text().splitlines()]


def validate_rows(rows):
    for row in rows:
        require(row["exit_code"] == 0 and row["valid"] and row["validity_checked"], "Invalid construction")
        require(row["vertices"] == row["expected_vertices"], "Vertex count mismatch")
        require(row["input_hash"] == row["output_hash"], "Point preservation hash mismatch")
        dimension = 1 if row["dataset"] in {"collinear", "exponential", "duplicates"} else 2
        require(row["dimension"] == dimension, "Unexpected effective dimension")
        for suffix in ("seconds", "cpu_seconds"):
            parts = [row[phase + "_" + suffix] for phase in ("sort", "insert", "total")]
            require(all(math.isfinite(value) and value >= 0 for value in parts), "Invalid time")
            require(math.isclose(parts[0] + parts[1], parts[2], rel_tol=1e-12, abs_tol=1e-12), "Phase sum mismatch")
        require(math.isfinite(row["insertion_path"]) and row["insertion_path"] > 0, "Invalid path")


def compare_csv(expected_path, actual_path):
    with expected_path.open(newline="") as stream:
        expected = list(csv.DictReader(stream))
    with actual_path.open(newline="") as stream:
        actual = list(csv.DictReader(stream))
    require(len(expected) == len(actual), "Summary row count mismatch")
    for a, b in zip(actual, expected):
        require(a.keys() == b.keys(), "Summary columns mismatch")
        for key in a:
            if a[key] != b[key]:
                try:
                    same = math.isclose(float(a[key]), float(b[key]), rel_tol=1e-13, abs_tol=1e-15)
                except ValueError:
                    same = False
                require(same, "Summary value mismatch: " + key)
    return actual


def main():
    schedule = json.loads((ROOT / "recorded_results/schedule.json").read_text())
    rows = samples("recorded_results")
    expected = Counter((case["dataset"], case["n"], case["seed"], case["variant"], trial)
                       for case in schedule["schedule"] for trial in range(case["repeats"]))
    actual = Counter((row["dataset"], row["n"], row["seed"], row["variant"], row["trial"]) for row in rows)
    require(actual == expected and all(value == 1 for value in actual.values()), "Incomplete or duplicate trials")
    require(len(rows) == 2730 and len(schedule["cases"]) == 46, "Campaign coverage mismatch")
    require(len(schedule["schedule"]) == 966 and schedule["seeds"] == list(range(1729, 1736)), "Schedule mismatch")
    validate_rows(rows)
    for directory in ("recorded_results", "validation_run"):
        instrument = json.loads((ROOT / directory / "instrumentation_validation.json").read_text())
        require(instrument["status"] == "PASS" and len(instrument["checks"]) == 32, "Instrumentation coverage mismatch")
        require(all(row["status"] == "PASS" and row["trials"] == 2 for row in instrument["checks"]), "Instrumentation failure")
    analyze.RESULTS = ROOT / "recorded_results"
    analyze.OUTPUT = ROOT / "rechecked_analysis"
    analyze.OUTPUT.mkdir(exist_ok=True)
    analyze.main()
    ratios = compare_csv(ROOT / "analysis/paired_ratios.csv", analyze.OUTPUT / "paired_ratios.csv")
    total = [row for row in ratios if row["numerator"] == "middle-corrected"
             and row["denominator"] == "middle-original" and row["metric"] == "total_seconds"]
    counts = {"cases": len(total), "point_estimate_below_one": sum(float(row["ratio"]) < 1 for row in total),
              "interval_entirely_above_one": sum(float(row["ci95_low"]) > 1 for row in total),
              "interval_entirely_below_one": sum(float(row["ci95_high"]) < 1 for row in total)}
    counts["interval_contains_one"] = counts["cases"] - counts["interval_entirely_above_one"] - counts["interval_entirely_below_one"]
    require(counts["cases"] == 46 and counts["point_estimate_below_one"] == 32
            and counts["interval_entirely_above_one"] == 7, "Issue summary mismatch")
    readme = (ROOT / "README.md").read_text()
    for row in total:
        line = (f'| {row["dataset"]} | {int(row["n"]):,} | {float(row["ratio"]):.4f} '
                f'| [{float(row["ci95_low"]):.4f}, {float(row["ci95_high"]):.4f}] |')
        require(line in readme, "README case missing or inconsistent: " + row["dataset"])
    focused = samples("validation_run")
    validate_rows(focused)
    lookup = {(row["seed"], row["variant"]): row for row in focused}
    require(len(focused) == 21 and len(lookup) == 21, "New York coverage mismatch")
    reference = {(row["seed"], row["variant"]): row for row in rows
                 if row["dataset"] == "roads_ny" and row["n"] == 264346}
    require(lookup.keys() == reference.keys(), "New York pairing mismatch")
    for key, row in lookup.items():
        old = reference[key]
        for field in ("input_hash", "input_order_hash", "insertion_order_hash", "vertices", "faces", "dimension"):
            require(row[field] == old[field], "New York geometry/hash mismatch")
        require(math.isclose(row["insertion_path"], old["insertion_path"], rel_tol=2e-14), "New York path mismatch")
    expected_order = [row for row in schedule["schedule"] if row["dataset"] == "roads_ny" and row["n"] == 264346]
    focused_environment = json.loads((ROOT / "validation_run/environment.json").read_text())
    require(focused_environment["schedule"] == expected_order, "New York relative schedule mismatch")
    with (ROOT / "validation_run/comparison.csv").open(newline="") as stream:
        comparison = list(csv.DictReader(stream))
    for row in comparison:
        metric = row["metric"]
        values = [lookup[seed, "middle-corrected"][metric] / lookup[seed, "middle-original"][metric]
                  for seed in schedule["seeds"]]
        computed = analyze.interval(values, "roads_ny-264346-middle-corrected-middle-original-" + metric)
        for key, value in zip(("corrected_over_original", "ci95_low", "ci95_high"), computed):
            require(math.isclose(value, float(row[key]), rel_tol=1e-13), "Focused New York statistic mismatch")
    result = {"status": "PASS", "original_constructions": len(rows), "policy_processes": 966,
              "paired_triplets": 322, "totals": counts, "focused_new_york_constructions": len(focused),
              "all_geometry_hashes_and_paths_match": True, "measurement_reruns": 0}
    (analyze.OUTPUT / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
