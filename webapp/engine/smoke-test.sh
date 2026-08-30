#!/usr/bin/env bash
#
# Smoke test for the compiled stroke_analysis binary.
#
# Verifies the C++ toolchain built correctly and that the binary produces its
# four expected output types from a CSV with the exact header the web app feeds
# it (SampleTimeFine,Acc_X,Acc_Y,Acc_Z,Gyr_X,Gyr_Y,Gyr_Z).
#
# Usage (inside the container):
#   docker compose exec web bash /app/engine/smoke-test.sh
#
# Usage (locally, after `make -C engine`):
#   bash engine/smoke-test.sh
#
# Optional: pass a real CSV to test with your own data instead of the synthetic one:
#   bash engine/smoke-test.sh /path/to/raw_xsens.csv   # header must already be Acc_X.. etc.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BINARY="${CPP_ENGINE_PATH:-$SCRIPT_DIR/stroke_analysis}"
WORK_DIR="$(mktemp -d)"
INPUT_CSV="$WORK_DIR/input.csv"
RESULT_DIR="$WORK_DIR/results"
CUSTOM_CSV="${1:-}"

pass() { printf '  \033[32mPASS\033[0m  %s\n' "$1"; }
fail() { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; FAILURES=$((FAILURES + 1)); }
info() { printf '\033[36m==>\033[0m %s\n' "$1"; }

FAILURES=0
cleanup() { rm -rf "$WORK_DIR"; }
trap cleanup EXIT

# --- 1. Binary present -------------------------------------------------------
info "Checking for compiled binary at: $BINARY"
if [[ ! -x "$BINARY" ]]; then
  fail "Binary not found or not executable. Build it first with: make -C engine"
  echo
  echo "The C++ sources may not be copied in yet. See engine/README.md."
  exit 1
fi
pass "Binary exists and is executable"

# --- 2. Prepare input CSV ----------------------------------------------------
mkdir -p "$RESULT_DIR"
if [[ -n "$CUSTOM_CSV" ]]; then
  info "Using provided CSV: $CUSTOM_CSV"
  if [[ ! -f "$CUSTOM_CSV" ]]; then
    fail "Provided CSV does not exist"
    exit 1
  fi
  cp "$CUSTOM_CSV" "$INPUT_CSV"
else
  info "Generating synthetic stroke CSV (10 strokes, 100 Hz)"
  # Synthesizes an Acc_Y sinusoid that crosses zero like real strokes so the
  # binary's peak/zero-crossing detector finds cycles.
  python3 - "$INPUT_CSV" <<'PY'
import math, sys
path = sys.argv[1]
fs, strokes, period = 100.0, 10, 1.2   # 100 Hz, 10 strokes, 1.2 s each
n = int(fs * strokes * period)
with open(path, "w") as f:
    f.write("SampleTimeFine,Acc_X,Acc_Y,Acc_Z,Gyr_X,Gyr_Y,Gyr_Z\n")
    for i in range(n):
        t = i / fs
        w = 2 * math.pi / period
        acc_y = 3.0 * math.sin(w * t)
        acc_x = 0.5 * math.sin(w * t + 0.4)
        acc_z = 9.81 + 0.3 * math.cos(w * t)
        gyr_x = 40.0 * math.cos(w * t)
        gyr_y = 15.0 * math.sin(w * t + 0.8)
        gyr_z = 10.0 * math.cos(w * t + 0.2)
        ts = int(t * 1e6)
        f.write(f"{ts},{acc_x:.5f},{acc_y:.5f},{acc_z:.5f},{gyr_x:.5f},{gyr_y:.5f},{gyr_z:.5f}\n")
print(f"wrote {n} rows")
PY
fi
[[ -s "$INPUT_CSV" ]] && pass "Input CSV ready ($(wc -l < "$INPUT_CSV") lines)" || { fail "Input CSV empty"; exit 1; }

# --- 3. Run the binary -------------------------------------------------------
info "Running: $BINARY -i input.csv -o results"
STDOUT_LOG="$WORK_DIR/stdout.log"
STDERR_LOG="$WORK_DIR/stderr.log"
MPLBACKEND="${MPLBACKEND:-Agg}" "$BINARY" -i "$INPUT_CSV" -o "$RESULT_DIR" \
  >"$STDOUT_LOG" 2>"$STDERR_LOG"
EXIT_CODE=$?

echo "  ---- stdout ----"
sed 's/^/  /' "$STDOUT_LOG" | head -20
if [[ -s "$STDERR_LOG" ]]; then
  echo "  ---- stderr ----"
  sed 's/^/  /' "$STDERR_LOG" | head -20
fi

if [[ $EXIT_CODE -eq 0 ]]; then
  pass "Binary exited 0"
else
  fail "Binary exited with code $EXIT_CODE"
fi

# --- 4. Verify outputs -------------------------------------------------------
info "Checking expected outputs in results/"

CYCLES_CSV="$RESULT_DIR/stroke_full_cycles.csv"
if [[ -s "$CYCLES_CSV" ]]; then
  ROWS=$(($(wc -l < "$CYCLES_CSV") - 1))
  pass "stroke_full_cycles.csv present ($ROWS strokes detected)"
  echo "  ---- cycles header ----"
  head -1 "$CYCLES_CSV" | sed 's/^/  /'
  [[ $ROWS -ge 1 ]] || fail "No stroke rows detected — check detection thresholds"
else
  fail "stroke_full_cycles.csv missing (core output — detection failed)"
fi

if [[ -s "$RESULT_DIR/stroke_phase_metrics.xlsx" ]]; then
  pass "stroke_phase_metrics.xlsx present ($(wc -c < "$RESULT_DIR/stroke_phase_metrics.xlsx") bytes) — xlnt OK"
else
  fail "stroke_phase_metrics.xlsx missing — check libxlnt linkage"
fi

if [[ -s "$RESULT_DIR/strokes_overview.png" ]]; then
  pass "strokes_overview.png present — matplotlib OK"
else
  fail "strokes_overview.png missing — check matplotlib/Agg backend"
fi

PLOT_COUNT=$(find "$RESULT_DIR/stroke_plots" -name 'stroke_*.png' 2>/dev/null | wc -l | tr -d ' ')
if [[ "$PLOT_COUNT" -ge 1 ]]; then
  pass "stroke_plots/ contains $PLOT_COUNT per-stroke plot(s)"
else
  fail "stroke_plots/ has no plots"
fi

# --- Summary -----------------------------------------------------------------
echo
if [[ $FAILURES -eq 0 ]]; then
  printf '\033[32mOK  Smoke test passed - the engine is wired up correctly.\033[0m\n'
  exit 0
else
  printf '\033[31mXX  Smoke test failed with %d issue(s). See messages above.\033[0m\n' "$FAILURES"
  exit 1
fi
