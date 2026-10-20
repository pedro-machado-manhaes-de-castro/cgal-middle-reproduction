# CGAL middle-policy recursion

Reproduce an early-stop ordering failure in CGAL's 2D `Hilbert_sort_middle_2` and measure the Delaunay tradeoffs of a [candidate correction](candidate.patch). The package contains the [original and corrected headers](headers), [benchmark](bench.cpp), and results for all 46 cases. The correction remains open for maintainer review.

**Western USA is the largest real input: 6,262,104 coordinates. Only New York's projected coordinates are bundled.** Other real inputs are downloaded and projected on demand, with source and coordinate SHA-256 checks. Every dataset can be measured independently.

## Setup

Requirements: Linux, little-endian hardware, GCC/libstdc++, Git, CMake >= 3.22, Ninja, Boost, GMP, MPFR, GNU time, Python and NumPy. From the package directory:

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

Keep this CGAL checkout clean. The commands select the supplied headers through separate include directories; applying the patch is unnecessary.

## Measure one dataset

Choose one command. Each measures the complete selected input with original middle, corrected middle and median, using seven seeds (21 policy processes and 21 timed constructions). It first builds the programs and checks the ordering diagnostic and insertion instrumentation with ASan/UBSan. Missing inputs are downloaded automatically; only the selected source is fetched.

| Dataset and source | Coordinates | Measurement command |
| --- | ---: | --- |
| **[Western USA](https://www.diag.uniroma1.it/~challenge9/data/USA-road-d/USA-road-d.W.co.gz)** | **6,262,104** | `python reproduce_dataset.py --dataset roads_w` |
| [California/Nevada](https://www.diag.uniroma1.it/~challenge9/data/USA-road-d/USA-road-d.CAL.co.gz) | 1,890,815 | `python reproduce_dataset.py --dataset roads_cal` |
| [GSHHG European crop](https://ftp.soest.hawaii.edu/gshhg/gshhg-bin-2.3.7.zip) | 1,826,772 | `python reproduce_dataset.py --dataset shorelines_europe` |
| [New York](https://www.diag.uniroma1.it/~challenge9/data/USA-road-d/USA-road-d.NY.co.gz) | 264,346 | `python reproduce_dataset.py --dataset roads_ny` |

For a smaller recorded case, select its size. The Western USA command below measures 65,536 points; it still needs the full coordinate file to reproduce the original seeded shuffle and prefix selection.

```sh
python reproduce_dataset.py --dataset roads_w --size 65536
```

Use `python reproduce_dataset.py --list` to list all 46 dataset/size choices. Without `--size`, the largest recorded size is used. Synthetic inputs need no downloads:

<details>
<summary>Commands for the eight synthetic families</summary>

| Family | Default n | Measurement command |
| --- | ---: | --- |
| Uniform | 4,194,304 | `python reproduce_dataset.py --dataset uniform` |
| Clustered | 4,194,304 | `python reproduce_dataset.py --dataset clustered` |
| Concentrated | 1,048,576 | `python reproduce_dataset.py --dataset concentrated` |
| Grid | 1,048,576 | `python reproduce_dataset.py --dataset grid` |
| Collinear | 16,384 | `python reproduce_dataset.py --dataset collinear` |
| Duplicates | 1,048,576 | `python reproduce_dataset.py --dataset duplicates` |
| Exponential | 16,384 | `python reproduce_dataset.py --dataset exponential` |
| Graded | 1,048,576 | `python reproduce_dataset.py --dataset graded` |

</details>

Outputs go to `dataset_results/<dataset>-<size>/`, for example `dataset_results/roads_w-6262104/`:

| File | Contents |
| --- | --- |
| `samples.jsonl` | Each policy/trial: `insertion_path`, `sort_seconds`, `insert_seconds`, `total_seconds`, CPU times, peak memory and validity checks |
| `comparison.csv` | Paired corrected/original ratios and individual 95% intervals |
| `validation.json` | Coverage, pairing and agreement with recorded geometry/order hashes |

Run commands sequentially because they share build directories. Repeating an interrupted command resumes completed processes; use a fresh extraction for an independent repetition. Each selection keeps its recorded seeds, repetitions and relative process order, omitting other cases. New timings remain separate from the supplied results. In the recorded libstdc++ environment, add `--require-recorded-order` to require matching geometry, insertion-order hashes and paths. The compatibility command `python reproduce_ny.py` runs full New York with output in `ny_results`.

## Inputs and full reproduction

To reconstruct coordinates without compiling or running Delaunay, choose one command:

| Dataset | Coordinate-only command |
| --- | --- |
| Western USA | `python rebuild_data.py --dataset roads_w` |
| California/Nevada | `python rebuild_data.py --dataset roads_cal` |
| GSHHG European crop | `python rebuild_data.py --dataset shorelines_europe` |
| New York | `python rebuild_data.py --dataset roads_ny` |

These write `rebuilt_data/<dataset>.xy`, with raw and projected hashes checked. The measurement commands reuse verified files; New York uses the bundled `data/roads_ny.xy`. [datasets.json](datasets.json) specifies sources, counts, projection centers, selections, hashes and attribution; [prepare_data.py](prepare_data.py) implements the conversion. Coordinates are headerless little-endian float64 x,y pairs in metres, using spherical Lambert azimuthal equal-area projection with radius 6,371,008.8 m. Duplicates are retained; road edges are unused. GSHHG selects full-resolution levels 1-5 within longitude [-25,45], latitude [30,72]. A projection-hash mismatch stops verification and may reflect a math-library difference.

`python rebuild_data.py` rebuilds all four inputs. `python reproduce.py` executes the entire 966-process, 2,730-construction campaign and writes new measurements to `results`; allow an hour or longer. Individual dataset commands above avoid this full run.

## Six-point ordering diagnostic

```sh
cmake -S witness -B build-witness -G Ninja -DCGAL_DIR="$CGAL_BASE" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_RELEASE="-O1 -UNDEBUG" -DCMAKE_CXX_COMPILER=g++
cmake --build build-witness --parallel 1
./build-witness/example_original
./build-witness/example_corrected
./build-witness/contract_witness_original
./build-witness/contract_witness_corrected
```

The [example](witness/example.cpp) uses `(t,t)` with `64*t = 0,5,9,7,11,64`. Original output is `0 5 9 7 11 64`; corrected output is `0 5 7 9 11 64`. The anchors fix `[0,1]^2`; five points occupy one child and its lower-left child. The original guard stops before later cuts separate them. All coordinates are exact binary fractions.

The [diagnostic](witness/contract_witness.cpp) checks all 720 permutations through both the default-limit functor and `hilbert_sort` with explicit middle: 714 failures per entry point originally, zero after correction. [Original](witness/contract_base.json) and [corrected](witness/contract_revised.json) results are supplied. CGAL's [reference](https://doc.cgal.org/6.1/Spatial_sorting/classCGAL_1_1Hilbert__sort__2.html) describes recursive midpoint subdivision; whether the early stop is intended remains a maintainer decision. The affected header is identical in the tested commit and main snapshot `cbc3b81d4607c6f3a9d1f3a47d8709a2be3e6bfa`.

## All 46 cases

```sh
python check_recorded.py
```

This checks all 2,730 original constructions and the separate 21-construction New York run, including coverage, pairing, validity, vertices, dimensions, hashes, phase sums and recomputed statistics. It checks every table entry below and writes `rechecked_analysis/verification.json`. It uses supplied measurements and performs no new Delaunay timing.

Corrected/original total-time estimates are below one in **32 of 46 cases**. Individual 95% intervals are wholly below one in 22 cases, wholly above one in seven, and contain one in 17. An interval containing one does not establish equivalence.

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

The [complete table](analysis/paired_ratios.csv) includes path, sorting, insertion, total and CPU ratios, with median as a third policy; [raw records](recorded_results/samples.jsonl) and the [schedule](recorded_results/schedule.json) are included. The original full-New-York total ratio was 1.051547 [1.020776, 1.095786]; the [separate run](validation_run/comparison.csv) gave 0.951908 [0.850160, 1.012797]. All 21 geometry/order hashes and paths matched, but the slowdown did not recur. Full production insertion sequences were not archived; timing results depend on the environment.

## Measurement and interpretation

The [protocol](DESIGN.txt) defines generators, sizes, sampling, repetitions and checks. Seeds are 1729-1735. Technical trials are reduced to their median per seed, followed by the geometric mean of seven paired ratios and 20,000-resample percentile-bootstrap intervals. Intervals are exploratory and unadjusted for multiplicity.

The benchmark measures the open Euclidean path of the actual insertion sequence, sorting, insertion and their contiguous total. Data preparation, warm-up, path measurement, validation and destruction are outside the timed interval; peak memory covers the whole process. Limits are one CPU, 8 GiB address space and 600 seconds per policy process. The 32 assertion-enabled ASan/UBSan configurations check the insertion trace, independently summed path and distinct-point count.

Recorded environment: Linux x86_64, GCC 15.2.0/libstdc++, Python 3.14.4, NumPy 2.5.3; see [environment.json](recorded_results/environment.json). Delaunay uses `Exact_predicates_inexact_constructions_kernel`, sequential sorting and `-O3 -DNDEBUG`; the diagnostic uses `Simple_cartesian<double>` and `-O1 -UNDEBUG`. CGAL adds `-frounding-math`. Phase timings show associations and do not isolate point-location, topology or cache effects.

## Licensing

Supplementary code is GPL version 3 or later; CGAL notices are retained in both headers. The [DIMACS page](https://www.diag.uniroma1.it/~challenge9/download.shtml) declares its data public domain. GSHHG 2.3.7, by Paul Wessel and Walter H. F. Smith, is LGPL version 3 or later. See [licenses](licenses) and [dataset attribution](datasets.json).
