# C++ stroke_analysis engine

The Next.js server spawns this binary for every uploaded run:

```
stroke_analysis -i <normalized.csv> -o <result_dir>
```

## Drop in your sources

Copy from your C++ repo into this folder:

```
engine/
  src/stroke_analysis.cpp        # already here (from your repo)
  utilities/                     # copy your whole utilities/ folder
    filters.cpp  filters.h
    metrics.cpp
    signal_features.h
    stroke_metrics.h
    math_utils.h
    cxxopts.hpp
  vendor/
    matplotlibcpp.h              # header-only matplotlib-cpp
```

From your C++ repo root:

```bash
cp -r utilities  /path/to/dashboard/engine/utilities
cp metrics.cpp filters.cpp /path/to/dashboard/engine/utilities/ 2>/dev/null || true
mkdir -p /path/to/dashboard/engine/vendor
curl -L -o /path/to/dashboard/engine/vendor/matplotlibcpp.h \
  https://raw.githubusercontent.com/lava/matplotlib-cpp/master/matplotlibcpp.h
```

`src/stroke_analysis.cpp` here is a copy. If you change the analysis in your
repo, re-copy it so the container rebuilds with your latest logic.

## Build locally

```bash
make -C engine          # -> engine/stroke_analysis
engine/stroke_analysis -i sample.csv -o /tmp/results
```

Requires `g++` (C++17), `python3-dev`, `python3-numpy`, `python3-matplotlib`
and `libxlnt-dev`.

## Input contract

The server never passes the raw Movella export. `lib/imu-csv.ts` strips the
metadata rows and emits exactly:

```
SampleTimeFine,Acc_X,Acc_Y,Acc_Z,Gyr_X,Gyr_Y,Gyr_Z
```

which is what `find_idx()` in `stroke_analysis.cpp` looks for.

## Output contract

Written into `-o <result_dir>` and harvested by `lib/cpp-engine.ts`:

| File | Stored as |
| --- | --- |
| `stroke_full_cycles.csv` | rows in `stroke_cycles` |
| `stroke_phase_metrics.xlsx` | downloadable artifact |
| `strokes_overview.png` | downloadable artifact |
| `stroke_plots/stroke_N.png` | per-stroke downloadable artifact |

Artifacts are copied to `ARTIFACT_DIR/<runId>/` (default `/app/data/runs`) and
served through `/api/runs/<id>/artifacts`.

## Toggle

Set `USE_CPP_ENGINE=false` to skip the binary; uploads still store samples.
If the binary is missing the upload succeeds and reports the engine error
instead of failing.

## Smoke test

After building, verify the toolchain end-to-end before relying on uploads:

```bash
# inside the running container
docker compose exec web bash /app/engine/smoke-test.sh

# or locally after `make -C engine`
bash engine/smoke-test.sh
```

It synthesizes a 10-stroke CSV with the exact `SampleTimeFine,Acc_X..Gyr_Z`
header, runs the binary, and checks all four outputs (cycles CSV, xlsx, overview
plot, per-stroke plots). Pass your own file to test real data:
`bash engine/smoke-test.sh /path/to/data.csv` (header must already be normalized).
Each PASS/FAIL line points at the likely culprit — e.g. a missing `.xlsx`
means libxlnt isn't linked, a missing `.png` means matplotlib/Agg failed.
