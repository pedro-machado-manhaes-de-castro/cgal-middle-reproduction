# Table: paired Delaunay measurements in analysis/paired_ratios.csv.
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import os
import platform
import random
import resource
import signal
import subprocess
import time


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)
VARIANTS = {"middle-original": ("original", "middle"),
            "middle-corrected": ("corrected", "middle"), "median": ("original", "median")}
SEEDS = list(range(1729, 1736))
MEMORY_LIMIT = 8 * 1024 ** 3
CPU = min(os.sched_getaffinity(0))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(2 ** 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix + ".part")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def constrained():
    resource.setrlimit(resource.RLIMIT_AS, (MEMORY_LIMIT, MEMORY_LIMIT))
    os.sched_setaffinity(0, {CPU})


def call(name, command, timeout=1200):
    logfile = RESULTS / (name + ".log")
    print(json.dumps({"starting": name}), flush=True)
    started = time.monotonic()
    with logfile.open("w") as stream:
        process = subprocess.run([str(x) for x in command], cwd=ROOT, stdout=stream,
                                 stderr=subprocess.STDOUT, timeout=timeout)
    with (RESULTS / "commands.jsonl").open("a") as stream:
        stream.write(json.dumps({"name": name, "command": list(map(str, command)),
                                 "seconds": time.monotonic() - started, "exit_code": process.returncode}) + "\n")
    if process.returncode:
        raise RuntimeError(f"{name} failed: {logfile}")


def hardware_info():
    processor = Path("/proc/cpuinfo").read_text()
    memory = Path("/proc/meminfo").read_text()
    fields = dict(line.split(":", 1) for line in processor.splitlines() if ":" in line)
    fields = {key.strip(): value.strip() for key, value in fields.items()}
    memory_fields = dict(line.split(":", 1) for line in memory.splitlines())
    return {"cpu_model": fields.get("model name", "unknown"),
            "logical_cpus": os.cpu_count(),
            "memory_kib": int(memory_fields["MemTotal"].split()[0]),
            "virtualized": "hypervisor" in fields.get("flags", "").split()}


def build():
    configuration = ROOT / "run_config.json"
    if not configuration.exists():
        configuration = ROOT / "experiment.json"
    frozen = json.loads(configuration.read_text())
    for name, digest in frozen["files"].items():
        if sha256(ROOT / name) != digest:
            raise RuntimeError("Frozen input changed: " + name)
    base = Path(os.environ.get("CGAL_BASE", frozen["cgal_base"]))
    if sha256(base / "Spatial_sorting/include/CGAL/Hilbert_sort_middle_2.h") != frozen["base_header_sha256"]:
        raise RuntimeError("Pinned CGAL base header mismatch")
    call("configure", ["cmake", "-S", ROOT, "-B", ROOT / "build", "-G", "Ninja",
                       "-DCMAKE_BUILD_TYPE=Release", "-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG",
                       "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON", "-DCMAKE_DISABLE_FIND_PACKAGE_TBB=TRUE",
                       "-DCGAL_DIR=" + str(base), "-DCMAKE_CXX_COMPILER=g++"])
    call("build", ["cmake", "--build", ROOT / "build", "--parallel", "1"])
    call("configure-checks", ["cmake", "-S", ROOT, "-B", ROOT / "build-checks", "-G", "Ninja",
                              "-DCMAKE_BUILD_TYPE=Debug", "-DCMAKE_CXX_FLAGS_DEBUG=-O1 -g -UNDEBUG -DBENCH_TRACE_INSERTION -fsanitize=address,undefined -fno-omit-frame-pointer",
                              "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON", "-DCMAKE_DISABLE_FIND_PACKAGE_TBB=TRUE",
                              "-DCGAL_DIR=" + str(base), "-DCMAKE_CXX_COMPILER=g++"])
    call("build-checks", ["cmake", "--build", ROOT / "build-checks", "--parallel", "1"])
    actual_headers = {}
    for variant in ["original", "corrected"]:
        dependency = subprocess.check_output(["ninja", "-C", str(ROOT / "build"), "-t", "deps",
                                              f"CMakeFiles/bench_{variant}.dir/bench.cpp.o"], text=True)
        (RESULTS / (variant + "-dependencies.txt")).write_text(dependency)
        headers = [line.strip() for line in dependency.splitlines()
                   if line.strip().endswith("/Hilbert_sort_middle_2.h")]
        expected = ROOT / "headers" / variant / "CGAL/Hilbert_sort_middle_2.h"
        if len(headers) != 1 or Path(headers[0]).resolve() != expected.resolve():
            raise RuntimeError("Wrong middle header in compiled binary")
        actual_headers[variant] = {"path": headers[0], "sha256": sha256(headers[0])}
    atomic_json(RESULTS / "environment.json", {
        "platform": platform.platform(), "hardware": hardware_info(), "cpu_affinity": CPU,
        "compiler": subprocess.check_output(["g++", "--version"], text=True),
        "python": platform.python_version(), "headers": actual_headers,
        "binaries": {v: sha256(ROOT / "build" / ("bench_" + v)) for v in ["original", "corrected"]},
        "base_commit": frozen["base_commit"], "address_space_limit_bytes": MEMORY_LIMIT})


def instrument_checks():
    checks = []
    changed = 0
    for family in ["uniform", "clustered", "concentrated", "grid", "collinear", "duplicates", "exponential", "graded"]:
        median_outputs = []
        for executable, policy in [("original", "middle"), ("corrected", "middle"),
                                    ("original", "median"), ("corrected", "median")]:
            tag = family + "-" + executable + "-" + policy
            dump = RESULTS / ("check-" + tag + ".csv")
            command = [ROOT / "build-checks" / ("bench_" + executable), policy,
                       "synthetic:" + family, "257", "1729", "2", "1", dump]
            process = subprocess.run(list(map(str, command)), capture_output=True, text=True, timeout=60)
            (RESULTS / ("check-" + tag + ".log")).write_text(process.stdout + process.stderr)
            if process.returncode:
                raise RuntimeError("Instrumented/sanitized construction failed: " + tag)
            rows = [json.loads(line) for line in process.stdout.splitlines()]
            with dump.open() as stream:
                points = [(float(row["x"]), float(row["y"])) for row in csv.DictReader(stream)]
            length = math.fsum(math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(points, points[1:]))
            if len(points) != 257 or len(rows) != 2:
                raise RuntimeError("Missing insertion trace or timing trial")
            for row in rows:
                if not math.isclose(length, row["insertion_path"], rel_tol=2e-14, abs_tol=1e-14):
                    raise RuntimeError("Path differs from independently summed insertion trace")
                if not math.isclose(row["sort_seconds"] + row["insert_seconds"], row["total_seconds"], rel_tol=1e-12):
                    raise RuntimeError("Phase timer sum mismatch")
                if not row["valid"] or row["vertices"] != len(set(points)):
                    raise RuntimeError("Vertex coverage or Delaunay validity failure")
                if row["input_hash"] != row["output_hash"]:
                    raise RuntimeError("Input changed during sorting")
                changed += not math.isclose(row["input_path"], row["insertion_path"], rel_tol=1e-8)
            if rows[0]["insertion_order_hash"] != rows[1]["insertion_order_hash"]:
                raise RuntimeError("Repeat changed deterministic insertion sequence")
            if policy == "median":
                median_outputs.append(rows[0]["insertion_order_hash"])
            checks.append({"case": tag, "path_reference": length, "trials": len(rows),
                           "vertices": rows[0]["vertices"], "status": "PASS"})
        if len(set(median_outputs)) != 1:
            raise RuntimeError("Median output differs between builds")
    if changed == 0:
        raise RuntimeError("Validation requires a case that changes the input path")
    atomic_json(RESULTS / "instrumentation_validation.json", {
        "status": "PASS", "checks": checks, "nontrivial_reordering_trials": changed,
        "method": "Insertion trace recorded inside the insertion loop in an assertion-enabled ASan/UBSan build; independent Python path sum and distinct point count; median control in both builds."})


def plan():
    datasets = json.loads((ROOT / "data/datasets.json").read_text())["datasets"]
    for dataset in datasets:
        if sha256(ROOT / dataset["file"]) != dataset["sha256"]:
            raise RuntimeError("Dataset changed before campaign")
    cases = []
    for family in ["uniform", "clustered", "concentrated", "grid", "collinear", "duplicates", "exponential", "graded"]:
        sizes = [4096, 16384] if family in ["collinear", "exponential"] else [16384, 65536, 262144, 1048576]
        if family in ["uniform", "clustered"]:
            sizes.append(4194304)
        for n in sizes:
            cases.append({"dataset": family, "source": "synthetic:" + family, "kind": "synthetic", "n": n})
    for dataset in datasets:
        sizes = sorted(set([n for n in [65536, 262144, 1048576, 4194304] if n < dataset["n"]] + [dataset["n"]]))
        for n in sizes:
            cases.append({"dataset": dataset["name"], "source": dataset["file"], "kind": "real", "n": n,
                          "dataset_sha256": dataset["sha256"], "full_dataset": n == dataset["n"]})
    chooser = random.Random(8297361)
    blocks = [(case, seed) for case in cases for seed in SEEDS]
    chooser.shuffle(blocks)
    schedule = []
    for case, seed in blocks:
        variants = list(VARIANTS)
        chooser.shuffle(variants)
        for variant in variants:
            repeats = 1 if case["dataset"] in ["collinear", "exponential"] else 5 if case["n"] <= 65536 else 3 if case["n"] <= 262144 else 1
            schedule.append({**case, "seed": seed, "variant": variant, "repeats": repeats})
    frozen = {"cases": cases, "seeds": SEEDS, "schedule": schedule,
              "policy_processes": len(schedule), "timed_trials": sum(x["repeats"] for x in schedule)}
    path = RESULTS / "schedule.json"
    if path.exists() and json.loads(path.read_text()) != frozen:
        raise RuntimeError("Attempted to change the frozen measurement schedule")
    atomic_json(path, frozen)
    return schedule


def execute(schedule):
    jobs = RESULTS / "jobs"
    jobs.mkdir(exist_ok=True)
    for index, case in enumerate(schedule):
        identifier = f"{case['dataset']}-{case['n']}-{case['seed']}-{case['variant']}"
        completed = jobs / (identifier + ".json")
        if completed.exists():
            continue
        atomic_json(RESULTS / "progress.json", {"state": "RUNNING", "index": index,
                                               "scheduled": len(schedule), "current": case})
        print(json.dumps({"index": index, "scheduled": len(schedule), "case": identifier}), flush=True)
        executable, policy = VARIANTS[case["variant"]]
        source = case["source"] if case["kind"] == "synthetic" else str(ROOT / case["source"])
        memory_file = jobs / (identifier + ".time")
        command = ["/usr/bin/time", "-f", "%M %U %S", "-o", str(memory_file),
                   str(ROOT / "build" / ("bench_" + executable)), policy, source,
                   str(case["n"]), str(case["seed"]), str(case["repeats"]), "1"]
        start = time.monotonic()
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                                   start_new_session=True, preexec_fn=constrained)
        try:
            stdout, stderr = process.communicate(timeout=600)
            code = process.returncode
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            stdout, stderr = process.communicate()
            code = 124
            stderr += " Per-policy 600-second process limit reached"
        rows = [json.loads(line) for line in stdout.splitlines() if line.startswith("{")]
        usage = memory_file.read_text().strip() if memory_file.exists() else None
        result = {"case": case, "exit_code": code, "elapsed_seconds": time.monotonic() - start,
                  "usage": usage, "samples": rows, "stderr": stderr, "command": command}
        if code == 0 and len(rows) != case["repeats"]:
            result["exit_code"] = 97
            result["stderr"] += " Missing timing trials"
        atomic_json(completed, result)
        if result["exit_code"]:
            print(json.dumps({"failure": identifier, "exit_code": result["exit_code"]}), flush=True)
    records = [json.loads(path.read_text()) for path in sorted(jobs.glob("*.json"))]
    with (RESULTS / "samples.jsonl").open("w") as stream:
        for record in records:
            for sample in record["samples"]:
                stream.write(json.dumps({**sample, **record["case"], "exit_code": record["exit_code"],
                                         "peak_rss_kib": int(record["usage"].split()[0]) if record["exit_code"] == 0 else None}) + "\n")
    failed = [record["case"] for record in records if record["exit_code"]]
    atomic_json(RESULTS / "progress.json", {"state": "COMPLETE_WITH_FAILURES" if failed else "COMPLETE",
                                           "policy_processes": len(records), "scheduled": len(schedule),
                                           "failed": failed})
    return failed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["prepare", "run", "all"], default="all")
    args = parser.parse_args()
    os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                      ASAN_OPTIONS="detect_leaks=1:halt_on_error=1", UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1")
    if args.phase in ["prepare", "all"]:
        build()
        instrument_checks()
        call("prepare-data", ["python", ROOT / "prepare_data.py"], timeout=3600)
    if args.phase in ["run", "all"]:
        if json.loads((RESULTS / "instrumentation_validation.json").read_text())["status"] != "PASS":
            raise RuntimeError("Instrumentation validation has not passed")
        failed = execute(plan())
        raise SystemExit(bool(failed))


if __name__ == "__main__":
    main()
