# Table: focused New York comparison in validation_run/comparison.csv.
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import os
import platform
import subprocess

import analyze
import run_campaign as campaign


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "ny_results"
BASE_COMMIT = "0490d81af4220f33005914676f362f349933fbc1"
POLICIES = ("middle-original", "middle-corrected", "median")
GEOMETRY = ("input_hash", "input_order_hash", "insertion_order_hash", "vertices", "faces", "dimension")


def checked(command):
    return subprocess.check_output(list(map(str, command)), text=True).strip()


def configure(source, build, flags, mode):
    campaign.call("configure-" + build, [
        "cmake", "-S", source, "-B", ROOT / build, "-G", "Ninja",
        "-DCMAKE_BUILD_TYPE=" + mode,
        "-DCMAKE_CXX_FLAGS_" + mode.upper() + "=" + flags,
        "-DCMAKE_DISABLE_FIND_PACKAGE_TBB=TRUE",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DCGAL_DIR=" + str(Path(os.environ["CGAL_BASE"]).resolve()),
        "-DCMAKE_CXX_COMPILER=g++",
    ])
    campaign.call("compile-" + build, ["cmake", "--build", ROOT / build, "--parallel", "1"])


def verify_sources():
    base = Path(os.environ["CGAL_BASE"]).resolve()
    if checked(["git", "-C", base, "rev-parse", "HEAD"]) != BASE_COMMIT:
        raise RuntimeError("Use the exact CGAL commit from README.md")
    subprocess.run(["git", "-C", str(base), "diff", "--exit-code", "HEAD", "--"], check=True)
    expected = json.loads((ROOT / "experiment.json").read_text())
    for name, digest in expected["files"].items():
        if campaign.sha256(ROOT / name) != digest:
            raise RuntimeError("Package source changed: " + name)
    source = (base / "Spatial_sorting/include/CGAL/Hilbert_sort_middle_2.h").read_bytes()
    original = (ROOT / "headers/original/CGAL/Hilbert_sort_middle_2.h").read_bytes()
    if source.replace(b"\r\n", b"\n") != original.replace(b"\r\n", b"\n"):
        raise RuntimeError("The original header differs from the pinned CGAL source")
    dataset = next(x for x in json.loads((ROOT / "datasets.json").read_text())["datasets"]
                   if x["name"] == "roads_ny")
    if campaign.sha256(ROOT / dataset["file"]) != dataset["sha256"]:
        raise RuntimeError("The New York projected coordinates have changed")
    return dataset


def verify_dependencies():
    for directory in ("build", "build-checks"):
        for variant in ("original", "corrected"):
            dependencies = checked([
                "ninja", "-C", ROOT / directory, "-t", "deps",
                f"CMakeFiles/bench_{variant}.dir/bench.cpp.o",
            ])
            headers = [line.strip() for line in dependencies.splitlines()
                       if line.strip().endswith("/Hilbert_sort_middle_2.h")]
            expected = ROOT / "headers" / variant / "CGAL/Hilbert_sort_middle_2.h"
            if len(headers) != 1 or Path(headers[0]).resolve() != expected.resolve():
                raise RuntimeError("The executable used the wrong header")


def check_witness():
    output = {}
    for variant in ("original", "corrected"):
        example = checked([ROOT / "build-witness" / ("example_" + variant)])
        diagnostic = json.loads(checked([ROOT / "build-witness" / ("contract_witness_" + variant)]))
        expected = json.loads((ROOT / "witness" / (
            "contract_base.json" if variant == "original" else "contract_revised.json"
        )).read_text())
        fields = ("example_output_numerators", "permutations", "direct_failures", "function_failures")
        if set(diagnostic) != set(fields) or diagnostic != {field: expected[field] for field in fields}:
            raise RuntimeError("The exhaustive diagnostic differs from the recorded result")
        sequence = [int(value) for value in example.split()]
        if sequence != diagnostic["example_output_numerators"]:
            raise RuntimeError("The issue listing differs from the exhaustive diagnostic")
        output[variant] = {"example": sequence, "diagnostic": diagnostic}
    campaign.atomic_json(RESULTS / "witness_validation.json", {"status": "PASS", **output})


def check_results(schedule, require_recorded_order):
    rows = [json.loads(line) for line in (RESULTS / "samples.jsonl").read_text().splitlines()]
    expected = {(x["seed"], x["variant"]) for x in schedule}
    if len(rows) != 21 or {(x["seed"], x["variant"]) for x in rows} != expected:
        raise RuntimeError("Expected all 21 full-New-York measurements")
    original = [json.loads(line) for line in (ROOT / "recorded_results/samples.jsonl").read_text().splitlines()]
    recorded = {(x["seed"], x["variant"]): x for x in original
                if x["dataset"] == "roads_ny" and x["n"] == 264346}
    lookup = {(x["seed"], x["variant"]): x for x in rows}
    order_matches = 0
    for row in rows:
        if row["exit_code"] or not row["validity_checked"] or not row["valid"]:
            raise RuntimeError("Failed triangulation validation")
        if row["vertices"] != row["expected_vertices"] or row["dimension"] != 2:
            raise RuntimeError("Unexpected vertex count or dimension")
        if row["input_hash"] != row["output_hash"]:
            raise RuntimeError("Sorting changed the point multiset")
        if not math.isclose(row["sort_seconds"] + row["insert_seconds"], row["total_seconds"], rel_tol=1e-12):
            raise RuntimeError("Wall time phases do not sum to the total")
        reference = recorded[row["seed"], row["variant"]]
        matches = all(row[field] == reference[field] for field in GEOMETRY)
        matches = matches and math.isclose(row["insertion_path"], reference["insertion_path"], rel_tol=2e-14)
        order_matches += matches
    for seed in range(1729, 1736):
        for field in ("input_hash", "input_order_hash", "vertices", "dimension"):
            if len({lookup[seed, policy][field] for policy in POLICIES}) != 1:
                raise RuntimeError("Unpaired inputs or triangulations")
    if require_recorded_order and order_matches != 21:
        raise RuntimeError("Insertion-order hashes, geometry or path differ from the recorded environment")
    comparisons = []
    for metric in ("insertion_path", "sort_seconds", "insert_seconds", "total_seconds"):
        values = [lookup[seed, "middle-corrected"][metric] / lookup[seed, "middle-original"][metric]
                  for seed in range(1729, 1736)]
        estimate, low, high = analyze.interval(
            values, f"roads_ny-264346-middle-corrected-middle-original-{metric}"
        )
        comparisons.append({"metric": metric, "corrected_over_original": estimate,
                            "ci95_low": low, "ci95_high": high})
    with (RESULTS / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparisons[0]))
        writer.writeheader()
        writer.writerows(comparisons)
    validation = {"status": "PASS", "timed_constructions": len(rows), "paired_seeds": 7,
                  "n": 264346, "recorded_geometry_matches": order_matches,
                  "timing_interpretation": "New measurements, separate from the original campaign; timing ratios need not repeat."}
    campaign.atomic_json(RESULTS / "validation.json", validation)
    print(json.dumps({**validation, "comparisons": comparisons}, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--require-recorded-order", action="store_true")
    args = parser.parse_args()
    os.chdir(ROOT)
    RESULTS.mkdir(exist_ok=True)
    campaign.RESULTS = RESULTS
    os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                      ASAN_OPTIONS="detect_leaks=1:halt_on_error=1",
                      UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1")
    dataset = verify_sources()
    configure(ROOT, "build", "-O3 -DNDEBUG", "Release")
    configure(ROOT, "build-checks", "-O1 -g -UNDEBUG -DBENCH_TRACE_INSERTION -fsanitize=address,undefined -fno-omit-frame-pointer", "Debug")
    configure(ROOT / "witness", "build-witness", "-O1 -UNDEBUG", "Release")
    verify_dependencies()
    check_witness()
    campaign.instrument_checks()
    historical = json.loads((ROOT / "recorded_results/schedule.json").read_text())
    schedule = [x for x in historical["schedule"] if x["dataset"] == "roads_ny" and x["n"] == 264346]
    identity = {"base_commit": BASE_COMMIT, "dataset_sha256": dataset["sha256"],
                "compiler": checked(["g++", "--version"]), "platform": platform.platform(),
                "binaries": {v: campaign.sha256(ROOT / "build" / ("bench_" + v))
                             for v in ("original", "corrected")},
                "schedule": schedule, "cpu_affinity": campaign.CPU}
    identity_path = RESULTS / "run_identity.json"
    if identity_path.exists() and json.loads(identity_path.read_text()) != identity:
        raise RuntimeError("Existing results belong to another build or schedule; use a fresh extraction")
    campaign.atomic_json(identity_path, identity)
    if campaign.execute(schedule):
        raise RuntimeError("The New York run contains failed processes")
    check_results(schedule, args.require_recorded_order)


if __name__ == "__main__":
    main()
