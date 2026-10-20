# CGAL middle-policy recursion

This repository reproduces a six-point ordering failure in `Hilbert_sort_middle_2` and the measured Delaunay tradeoffs of a candidate correction. It targets CGAL commit `0490d81af4220f33005914676f362f349933fbc1`. The candidate affects middle ordering in 2D. Its suitability for CGAL remains open for review.

## Evidence and commands

| Question | Evidence | Reproduction |
| --- | --- | --- |
| What stops refinement, and what changes? | [Original header](headers/original/CGAL/Hilbert_sort_middle_2.h), [candidate patch](candidate.patch), [candidate header](headers/corrected/CGAL/Hilbert_sort_middle_2.h) | Build the six-point diagnostic below |
| Does the six-point example fail across input orders? | [Example](witness/example.cpp), [exhaustive diagnostic](witness/contract_witness.cpp), [original result](witness/contract_base.json), [candidate result](witness/contract_revised.json) | Build the diagnostic, or `python reproduce_ny.py` |
| Which 46 cases support the performance statements? | [All cases below](#all-46-cases), [paired ratios](analysis/paired_ratios.csv), [all 2,730 measurements](recorded_results/samples.jsonl), [schedule](recorded_results/schedule.json) | `python check_recorded.py` |
| Where do all real and synthetic inputs come from? | [Public sources and projection hashes](datasets.json), [projection code](prepare_data.py), [synthetic generators](bench.cpp) | `python rebuild_data.py`; `python reproduce.py` generates synthetic inputs |
| How are paths, times and validity measured? | [Benchmark](bench.cpp), [measurement protocol](DESIGN.txt), [instrumentation checks](recorded_results/instrumentation_validation.json), [environment](recorded_results/environment.json) | `python reproduce_ny.py`; `python reproduce.py` for all cases |
| Was the New York slowdown repeated? | [Separate observations](validation_run/samples.jsonl), [ratios](validation_run/comparison.csv), [environment and order](validation_run/environment.json) | `python check_recorded.py` checks these observations; `python reproduce_ny.py` makes new measurements |

## Setup

Use Linux, little-endian hardware, GCC/libstdc++, Git, CMake 3.22 or later, Ninja, Boost, GMP, MPFR, GNU time, Python and NumPy. The recorded environment used GCC 15.2.0, Python 3.14.4 and NumPy 2.5.3. The benchmark uses `Exact_predicates_inexact_constructions_kernel`; the six-point diagnostic uses `Simple_cartesian<double>`. Sorting is sequential. Delaunay measurements use `-O3 -DNDEBUG`; the diagnostic uses `-O1 -UNDEBUG`. CGAL adds `-frounding-math`.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install numpy
git init cgal
git -C cgal remote add origin https://github.com/CGAL/cgal.git
git -C cgal fetch --depth 1 origin 0490d81af4220f33005914676f362f349933fbc1
git -C cgal checkout --detach FETCH_HEAD
export CGAL_BASE="$(realpath cgal)"
sha256sum --check SHA256SUMS
```

Keep the baseline checkout clean. Separate include directories select the two supplied headers, so applying the patch is unnecessary for these commands. Sources and inputs are checked before compilation. `SHA256SUMS` identifies the distributed files. The benchmark body and both sorting headers are those used for the recorded measurements; the benchmark's opening comment names the table provided here.

## Six-point ordering diagnostic

```sh
cmake -S witness -B build-witness -G Ninja -DCGAL_DIR="$CGAL_BASE" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_RELEASE="-O1 -UNDEBUG" -DCMAKE_CXX_COMPILER=g++
cmake --build build-witness --parallel 1
./build-witness/example_original
./build-witness/example_corrected
./build-witness/contract_witness_original
./build-witness/contract_witness_corrected
```

The points are `(t,t)` with integer numerators `64*t = 0,5,9,7,11,64`. The original prints `0 5 9 7 11 64`; the candidate prints `0 5 7 9 11 64`. All coordinates are exact binary fractions. The anchors fix `[0,1]^2`; five points occupy `[0,1/2]^2` and then its lower-left child. The original recursion guard stops before later cuts can separate them.

| Entry point | Input permutations | Original failures | Candidate failures |
| --- | ---: | ---: | ---: |
| Direct functor, default leaf limit 1 | 720 | 714 | 0 |
| `hilbert_sort`, explicit middle | 720 | 714 | 0 |

The [CGAL 6.1 reference](https://doc.cgal.org/6.1/Spatial_sorting/classCGAL_1_1Hilbert__sort__2.html) and [manual](https://doc.cgal.org/6.1/Spatial_sorting/index.html) describe recursive midpoint subdivision. Whether the early stop is intended requires a maintainer decision. The affected header is identical at the tested commit and main snapshot `cbc3b81d4607c6f3a9d1f3a47d8709a2be3e6bfa` (Git blob `9c48ecbd8356bf391e0aa8e6f91847c5ae2cc3d8`).

## All 46 cases

`python check_recorded.py` checks every original trial and all seven paired seeds, recomputes the published ratios and intervals, checks this table, and independently recomputes the separate New York comparison. It checks dimensions, validity flags, vertex counts, point-preservation hashes, phase sums, exact trial coverage and the recorded insertion-trace checks. Generated tables and the verification result go to `rechecked_analysis`.

The corrected/original total-time point estimate is below one in **32 of 46 cases**. Individual intervals lie entirely below one in 22 cases, entirely above one in seven, and contain one in 17. Ratios are geometric means of seven paired seed ratios, after taking the median of technical trials within each seed. Individual percentile-bootstrap intervals use 20,000 resamples and a deterministic seed derived from the comparison label. They are exploratory and unadjusted for multiplicity. Containing one does not establish equivalence.

<details>
<summary>All corrected/original total-time ratios, including favorable and unfavorable cases</summary>

| Dataset | n | Corrected/original total | Individual 95% interval |
| --- | ---: | ---: | --- |
| uniform | 16,384 | 1.0081 | [1.0048, 1.0116] |
| uniform | 65,536 | 0.9990 | [0.9864, 1.0089] |
| uniform | 262,144 | 0.9957 | [0.9764, 1.0117] |
| uniform | 1,048,576 | 1.0138 | [0.9986, 1.0348] |
| uniform | 4,194,304 | 1.0426 | [0.9960, 1.1353] |
| clustered | 16,384 | 0.4973 | [0.4912, 0.5031] |
| clustered | 65,536 | 0.2464 | [0.2455, 0.2473] |
| clustered | 262,144 | 0.1134 | [0.1032, 0.1293] |
| clustered | 1,048,576 | 0.0520 | [0.0495, 0.0541] |
| clustered | 4,194,304 | 0.0104 | [0.0102, 0.0107] |
| concentrated | 16,384 | 0.3162 | [0.2078, 0.5872] |
| concentrated | 65,536 | 0.1006 | [0.0778, 0.1670] |
| concentrated | 262,144 | 0.1044 | [0.0405, 0.3110] |
| concentrated | 1,048,576 | 0.0075 | [0.0046, 0.0192] |
| grid | 16,384 | 0.9911 | [0.9752, 1.0005] |
| grid | 65,536 | 0.9960 | [0.9919, 0.9992] |
| grid | 262,144 | 0.9952 | [0.9759, 1.0158] |
| grid | 1,048,576 | 0.9975 | [0.9911, 1.0029] |
| collinear | 4,096 | 0.9721 | [0.9606, 0.9847] |
| collinear | 16,384 | 0.9984 | [0.9838, 1.0117] |
| duplicates | 16,384 | 0.8894 | [0.6886, 1.0167] |
| duplicates | 65,536 | 1.0180 | [1.0134, 1.0222] |
| duplicates | 262,144 | 1.0149 | [1.0095, 1.0204] |
| duplicates | 1,048,576 | 1.0056 | [0.9914, 1.0216] |
| exponential | 4,096 | 0.7043 | [0.6197, 0.7858] |
| exponential | 16,384 | 0.6723 | [0.6314, 0.7156] |
| graded | 16,384 | 1.0041 | [0.9846, 1.0184] |
| graded | 65,536 | 0.9965 | [0.9762, 1.0106] |
| graded | 262,144 | 1.0654 | [1.0038, 1.1707] |
| graded | 1,048,576 | 1.0809 | [0.9579, 1.2409] |
| roads_ny | 65,536 | 1.0159 | [1.0075, 1.0282] |
| roads_ny | 262,144 | 1.0140 | [0.9993, 1.0300] |
| roads_ny | 264,346 | 1.0515 | [1.0208, 1.0958] |
| roads_cal | 65,536 | 1.0062 | [1.0004, 1.0125] |
| roads_cal | 262,144 | 1.0037 | [0.9833, 1.0264] |
| roads_cal | 1,048,576 | 0.9873 | [0.9685, 1.0010] |
| roads_cal | 1,890,815 | 0.9833 | [0.9620, 1.0012] |
| roads_w | 65,536 | 0.9358 | [0.9301, 0.9402] |
| roads_w | 262,144 | 0.8091 | [0.7451, 0.8487] |
| roads_w | 1,048,576 | 0.7194 | [0.7124, 0.7275] |
| roads_w | 4,194,304 | 0.4987 | [0.4850, 0.5079] |
| roads_w | 6,262,104 | 0.4123 | [0.3940, 0.4373] |
| shorelines_europe | 65,536 | 0.9906 | [0.9840, 0.9966] |
| shorelines_europe | 262,144 | 0.9700 | [0.9611, 0.9789] |
| shorelines_europe | 1,048,576 | 0.9253 | [0.9158, 0.9338] |
| shorelines_europe | 1,826,772 | 0.8815 | [0.8693, 0.8909] |

</details>

The [complete CSV](analysis/paired_ratios.csv) also contains sorting, insertion, path and CPU ratios, with median as a third policy. Filter `numerator=middle-corrected` and `denominator=middle-original` for the candidate comparison. For example, `concentrated`, n=1,048,576 has sorting ratio 3.0100 and contiguous total ratio 0.0075. Phase timing supports associations; it does not separate the causal contributions of point location, topology and caching.

## Inputs and full reproduction

| Dataset | Full coordinate count | Public source | Projection center (longitude, latitude) |
| --- | ---: | --- | --- |
| `roads_ny` | 264,346 | [DIMACS New York](https://www.diag.uniroma1.it/~challenge9/data/USA-road-d/USA-road-d.NY.co.gz) | (-74.0, 40.8) |
| `roads_cal` | 1,890,815 | [DIMACS California/Nevada](https://www.diag.uniroma1.it/~challenge9/data/USA-road-d/USA-road-d.CAL.co.gz) | (-119.5, 37.25) |
| `roads_w` | 6,262,104 | [DIMACS Western USA](https://www.diag.uniroma1.it/~challenge9/data/USA-road-d/USA-road-d.W.co.gz) | (-115.0, 38.5) |
| `shorelines_europe` | 1,826,772 | [GSHHG 2.3.7](https://ftp.soest.hawaii.edu/gshhg/gshhg-bin-2.3.7.zip), full-resolution vertices in [-25,45] longitude and [30,72] latitude, levels 1-5 | (10.0, 52.0) |

All use spherical Lambert azimuthal equal-area projection, radius 6,371,008.8 metres. Repeated coordinates and original record order are retained before seeded shuffling; graph edges are unused. The projected format is headerless little-endian float64 x,y pairs. Exact New York coordinates are included in `data/roads_ny.xy` (4,229,536 bytes). The other inputs are regenerated from their public sources to keep the repository small. [datasets.json](datasets.json) gives the raw and projected SHA-256 hashes, counts, selections and attribution for all four.

```sh
python rebuild_data.py
```

This independently downloads and regenerates all four datasets in `rebuilt_data`, checks raw and projected hashes, and preserves the included input. To check one source, add `--dataset roads_w` (or another table identifier). Projection hashes can depend on math-library versions; a mismatch stops verification. Eight synthetic families, `uniform`, `clustered`, `concentrated`, `grid`, `collinear`, `duplicates`, `exponential` and `graded`, are defined exactly in `bench.cpp`; no synthetic data files are needed.

To make new measurements for every case, use a fresh extraction and run:

```sh
python reproduce.py
```

This builds both variants, checks the included headers, performs 32 ASan/UBSan instrumentation checks (64 constructions), regenerates the real inputs with their recorded hashes, and executes the same 966-process schedule with 2,730 timed constructions. It uses seeds 1729 through 1735, up to 6,262,104 coordinates, 8 GiB address space per process, one CPU and a 600-second process limit. It may take an hour or longer. [DESIGN.txt](DESIGN.txt) specifies every size, repetition count, warm-up, input ordering and time boundary. An interrupted run resumes its completed processes. Fresh measurements go to `results`; a fresh extraction is required for an independent run.

## Focused New York measurements

```sh
python reproduce_ny.py
```

This builds and runs the ordering diagnostic, the 32 ASan/UBSan configurations, and all 264,346 New York coordinates for three policies and seven seeds (21 constructions). It preserves the relative New York process order from the full schedule, omitting intervening cases. It writes `ny_results/samples.jsonl`, `comparison.csv` and `validation.json`. The fields `insertion_path`, `sort_seconds`, `insert_seconds` and `total_seconds` give the open Euclidean path and the three wall-time intervals. Corresponding CPU times and peak process memory are also recorded. Preparation, warm-up, path measurement, validity checks and destruction lie outside the timing interval.

The original full-campaign total-time ratio was 1.051547 [1.020776, 1.095786]. The separate 21-construction run in [validation_run](validation_run) gave 0.951908 [0.850160, 1.012797]. Input and insertion-order hashes, vertex/face counts, dimensions and path lengths matched all 21 original records. Full insertion sequences were not archived for these production measurements. The second run did not reproduce the slowdown, and the two runs are kept separate. New timing results can differ in magnitude or sign. With the recorded libstdc++ environment, add `--require-recorded-order` to enforce agreement of recorded geometry, order hashes and paths.

The ordering diagnostic and New York procedure were tested against a clean official checkout. The complete original campaign and its separate New York rerun are supplied. Checking all recorded statistics does not constitute a new execution of all 966 policy processes.

## Licensing

CGAL notices are retained in both headers. Supplementary code is GPL version 3 or later; see [licenses](licenses). The [DIMACS download page](https://www.diag.uniroma1.it/~challenge9/download.shtml) declares its data public domain. GSHHG 2.3.7 is by Paul Wessel and Walter H. F. Smith under LGPL version 3 or later. Attribution and data licenses are recorded in `datasets.json`.
