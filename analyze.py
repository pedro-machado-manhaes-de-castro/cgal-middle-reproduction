# Table: all paired Delaunay comparisons in analysis/paired_ratios.csv.
from collections import defaultdict
from pathlib import Path
import csv
import hashlib
import json
import math
import statistics

import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUTPUT = ROOT / "analysis"
OUTPUT.mkdir(exist_ok=True)
POLICIES = ["middle-original", "middle-corrected", "median"]
METRICS = ["sort_seconds", "insert_seconds", "total_seconds", "insertion_path",
           "sort_cpu_seconds", "insert_cpu_seconds", "total_cpu_seconds"]


def write_csv(name, rows):
    with (OUTPUT / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def interval(values, label):
    values = np.asarray(values, dtype=float)
    logs = np.log(values)
    seed = int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], "little")
    random = np.random.default_rng(seed)
    draws = random.integers(0, len(values), size=(20000, len(values)))
    means = np.exp(logs[draws].mean(axis=1))
    low, high = np.quantile(means, [0.025, 0.975])
    return float(np.exp(logs.mean())), float(low), float(high)


def main():
    progress = json.loads((RESULTS / "progress.json").read_text())
    schedule = json.loads((RESULTS / "schedule.json").read_text())
    rows = [json.loads(line) for line in (RESULTS / "samples.jsonl").read_text().splitlines()]
    expected = {(x["dataset"], x["n"], x["seed"], x["variant"]): x["repeats"] for x in schedule["schedule"]}
    groups = defaultdict(list)
    for row in rows:
        key = (row["dataset"], row["n"], row["seed"], row["variant"])
        groups[key].append(row)
        if row["input_hash"] != row["output_hash"] or not row["valid"] or not row["validity_checked"]:
            raise RuntimeError("Failed preservation or validity check")
        if row["vertices"] != row["expected_vertices"]:
            raise RuntimeError("Distinct-point coverage mismatch")
        if not math.isclose(row["sort_seconds"] + row["insert_seconds"], row["total_seconds"], rel_tol=1e-12):
            raise RuntimeError("Phase times do not sum to total")
    incomplete = [key for key, count in expected.items() if len(groups[key]) != count]
    if progress["failed"] or incomplete:
        (OUTPUT / "INCOMPLETE.json").write_text(json.dumps({"failed": progress["failed"], "incomplete": incomplete}, indent=2))
        raise RuntimeError("Campaign contains failures or incomplete pairs; preserve censored cases")
    seeds = []
    for key, trials in sorted(groups.items()):
        for field in ["input_hash", "input_order_hash", "insertion_order_hash", "vertices", "dimension", "insertion_path"]:
            if len({row[field] for row in trials}) != 1:
                raise RuntimeError("Technical trials changed geometry or insertion order")
        dataset, n, seed, variant = key
        seeds.append({"dataset": dataset, "kind": trials[0]["kind"], "n": n, "seed": seed, "variant": variant,
                      "trials": len(trials), **{metric: statistics.median(x[metric] for x in trials) for metric in METRICS},
                      "input_path": trials[0]["input_path"], "vertices": trials[0]["vertices"],
                      "input_hash": trials[0]["input_hash"], "input_order_hash": trials[0]["input_order_hash"],
                      "insertion_order_hash": trials[0]["insertion_order_hash"],
                      "peak_rss_kib": max(x["peak_rss_kib"] for x in trials)})
    write_csv("seed_measurements.csv", seeds)
    lookup = {(x["dataset"], x["n"], x["seed"], x["variant"]): x for x in seeds}
    for case in schedule["cases"]:
        for seed in schedule["seeds"]:
            triplet = [lookup[(case["dataset"], case["n"], seed, variant)] for variant in POLICIES]
            for field in ["input_hash", "input_order_hash", "vertices"]:
                if len({x[field] for x in triplet}) != 1:
                    raise RuntimeError("Policies did not use identical inputs and vertex sets")
    summaries = []
    ratios = []
    for case in schedule["cases"]:
        dataset, n = case["dataset"], case["n"]
        for variant in POLICIES:
            subset = [lookup[(dataset, n, seed, variant)] for seed in schedule["seeds"]]
            summaries.append({"dataset": dataset, "kind": case["kind"], "n": n, "variant": variant,
                              "seeds": len(subset), "vertices_min": min(x["vertices"] for x in subset),
                              "vertices_max": max(x["vertices"] for x in subset),
                              **{metric: statistics.median(x[metric] for x in subset) for metric in METRICS},
                              "peak_rss_kib": max(x["peak_rss_kib"] for x in subset)})
        for numerator, denominator in [("middle-corrected", "middle-original"), ("median", "middle-original"),
                                        ("middle-corrected", "median")]:
            for metric in METRICS:
                values = [lookup[(dataset, n, seed, numerator)][metric] / lookup[(dataset, n, seed, denominator)][metric]
                          for seed in schedule["seeds"]]
                estimate, low, high = interval(values, f"{dataset}-{n}-{numerator}-{denominator}-{metric}")
                ratios.append({"dataset": dataset, "kind": case["kind"], "n": n, "numerator": numerator,
                               "denominator": denominator, "metric": metric, "ratio": estimate,
                               "ci95_low": low, "ci95_high": high, "pairs": len(values),
                               "min_pair_ratio": min(values), "max_pair_ratio": max(values),
                               "lower_in_pairs": sum(x < 1 - 1e-12 for x in values)})
    write_csv("timings_and_paths.csv", summaries)
    write_csv("paired_ratios.csv", ratios)
    print(json.dumps({"status": "PASS", "cases": len(schedule["cases"]),
                      "policy_processes": len(groups), "timed_constructions": len(rows)}, indent=2))


if __name__ == "__main__":
    main()
