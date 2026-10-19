# CGAL middle-policy recursion: reproducing the ordering and timing results

This package supports a discussion of the stopping condition in Hilbert_sort_middle_2. It contains the six-point counterexample, the proposed header change, the exact New York input, and the measurements behind the reported performance tradeoffs. The candidate changes only middle ordering in 2D. Its release suitability remains a question for review.

## 1. Reproduce New York

Use Linux on a little-endian machine with GCC, Git, CMake (at least 3.22), Ninja, Boost, GMP, MPFR, GNU time, Python, NumPy and Matplotlib. The recorded measurements used GCC 15.2.0, libstdc++, Python 3.14.4, NumPy 2.5.3, and Matplotlib 3.11.2. A Python environment can be prepared with:

    python3 -m venv .venv
    . .venv/bin/activate
    python -m pip install numpy matplotlib

Obtain the exact CGAL revision in a separate, clean directory. These commands retrieve one commit from the official repository:

    git init cgal
    git -C cgal remote add origin https://github.com/CGAL/cgal.git
    git -C cgal fetch --depth 1 origin 0490d81af4220f33005914676f362f349933fbc1
    git -C cgal checkout --detach FETCH_HEAD
    export CGAL_BASE="$(realpath cgal)"

From the extracted package directory, run:

    sha256sum --check SHA256SUMS
    python reproduce_ny.py

The script builds the original and candidate versions using different include directories. It verifies which header each executable actually includes. The supplied candidate.patch is available for inspection; applying it to CGAL_BASE is unnecessary for this comparison, which requires a clean baseline checkout.

The script first compiles and checks the six-point example and all 720 permutations. It then checks the insertion instrumentation in 32 assertion-enabled ASan/UBSan configurations, covering all eight synthetic families, both builds and both policies. Finally it measures original middle, corrected middle and median on all 264,346 New York coordinates, using the original seven seeds (1729 through 1735), one timed construction per policy and seed, and the relative execution order from the recorded campaign. This gives 21 timed constructions and seven paired triplets. Each process has its original warm-up, CPU affinity, 8 GiB address-space limit and 600-second timeout. The C++ benchmark and both supplied headers are preserved byte for byte from the recorded campaign.

Results are written to ny_results/samples.jsonl, ny_results/comparison.csv and ny_results/validation.json. They include the actual open insertion path, sorting, insertion, contiguous total, CPU time and peak process memory. Input preparation, warm-up and validity checks are outside the timed interval. A second invocation resumes this same run and does not create new statistical repetitions. Use a fresh extraction for another timing run.

The six-point defect is a functional observation. A New York timing ratio is an observation in a particular environment. A new execution can have a different ratio or even a different sign of the timing difference. Standard-library implementations can also produce different permutations. With the recorded GCC/libstdc++ environment, add --require-recorded-order to require all 21 geometry and insertion-order comparisons to match the historical records. Independently of that option, every new run checks paired input permutations, Delaunay validity, distinct-vertex counts, dimension, point preservation and phase-time sums.

## 2. Identify and regenerate the New York input

The source is the New York City coordinate file from the 9th DIMACS Implementation Challenge:

    https://www.diag.uniroma1.it/~challenge9/download.shtml
    https://www.diag.uniroma1.it/~challenge9/data/USA-road-d/USA-road-d.NY.co.gz

All 264,346 coordinate records are used; graph edges are not used. The original integer longitude and latitude values are multiplied by 1e-6. The spherical Lambert azimuthal equal-area projection has center (-74.0, 40.8) degrees and radius 6,371,008.8 metres. Duplicate locations and the source record order are retained before the seeded shuffle. The benchmark builds planar Delaunay triangulations and measures Euclidean lengths in this projected plane.

The exact projected input is included as data/roads_ny.xy: little-endian float64 x,y pairs, no header, 4,229,536 bytes. Its SHA-256 is:

    4c9bc16a9176cf56aadcd0536aa759c5e5c151b92c1175bf0cb5e8e0ea43c571

The compressed DIMACS source has 2,058,598 bytes and SHA-256:

    02547a628742164e02bd73f59e33a65d23dfac8fc54c09275c8b9169c38eafdf

To download, verify and independently regenerate the projection:

    python rebuild_ny.py

This writes to rebuilt_ny and leaves the benchmark input intact. Projection hashes can vary with floating-point math libraries; the included projected file fixes the exact input used for the reported measurements. The projection routine is in prepare_data.py, with parameters and source provenance in datasets.json.

## 3. Check the reported numbers without rerunning triangulations

recorded_results/samples.jsonl contains all 2730 original timed constructions, including every dataset and size, all three policies and all seven seeds. The full randomized schedule is in recorded_results/schedule.json. Environment information, instrumentation checks and the data dictionary accompany these records. They support all 46 dataset/size comparisons cited in the issue, including Western USA and uniform points.

Run:

    python check_recorded.py

This rechecks coverage, pairing, validity, vertex counts and phase times, regenerates the four numerical CSV tables, and compares them with analysis. The unit of analysis is a seed, using the median of its technical repetitions. Ratios are geometric means of seven paired ratios. Individual 95% percentile-bootstrap intervals use 20,000 resamples and a deterministic seed derived from the comparison label. These exploratory intervals have no multiplicity adjustment and do not establish equivalence when they contain one.

For full New York, the original corrected/original total-time ratio was 1.051547 with interval [1.020776, 1.095786]. The original path ratio was 0.995036. These numbers can be recovered from analysis/paired_ratios.csv by selecting roads_ny, n=264346, middle-corrected over middle-original. Each original per-seed row is also available in analysis/seed_measurements.csv. The newly measured comparison remains separate in ny_results.

The New York procedure was tested on 19 October 2026 from a fresh extraction, with a clean checkout fetched from the official CGAL repository and GCC 15.2.0. All 21 New York insertion sequences and paths matched the original records. The new total-time ratio was 0.951908 [0.850160, 1.012797], which did not reproduce the earlier slowdown. This focused run retained the relative order of the 21 New York processes and omitted the other cases interleaved in the original campaign. Its raw measurements and checks are included separately in validation_run. It validates the reproduction procedure and illustrates timing variability; the two campaigns are not pooled as independent observations.

To repeat the complete 46-case measurement campaign, use a fresh extraction, set CGAL_BASE, and run python reproduce.py. That command downloads the remaining public datasets, verifies their recorded projected hashes, and runs the frozen 966-process schedule. It can take about an hour or longer and requires more disk space and memory than the focused New York run. The full rerun was performed for the original study; the validation of this distribution repeats New York and reanalyses all recorded results.

## 4. Reproduce only the ordering defect

The issue listing is witness/example.cpp. witness/contract_witness.cpp checks every permutation through both the direct functor and hilbert_sort. These can be built separately:

    cmake -S witness -B build-witness -G Ninja -DCGAL_DIR="$CGAL_BASE" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_RELEASE="-O1 -UNDEBUG" -DCMAKE_CXX_COMPILER=g++
    cmake --build build-witness --parallel 1
    ./build-witness/example_original
    ./build-witness/example_corrected
    ./build-witness/contract_witness_original
    ./build-witness/contract_witness_corrected

The example prints 0 5 9 7 11 64 with the original header and 0 5 7 9 11 64 with the candidate. The exhaustive diagnostic records 714 failures out of 720 permutations in each original entry point and zero with the candidate. Its historical JSON outputs are included beside the sources.

## Source notices and data permissions

Original CGAL notices remain in both supplied headers. The supplementary benchmark and reproduction code is provided under GNU GPL version 3 or later; see licenses/COPYINGv3. DIMACS states that the coordinate files are public domain on its download page. The optional full campaign also uses GSHHG 2.3.7 under LGPL version 3 or later; its notices and license are in licenses. datasets.json records the attribution and sources. CGAL and all other dependencies retain their respective licenses.
