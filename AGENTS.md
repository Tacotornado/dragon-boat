# AGENTS.md — Dragon Boat / Canoe IMU Analytics Project

This file is read automatically by Antigravity before every task in this workspace.
It captures durable project context, decisions already made, and conventions —
so you don't need to re-explain the project in every new conversation.
Task-specific specs (e.g. "implement module X") are separate files; this file is the
standing background those specs build on.

---

## 1. What this project is

A real-time performance-analytics system for dragon boat / canoe racers, built around
Xsens DOT IMU sensors. Entry for the 2026 "Zhongji Technology Application Creative
Competition" (theme: Sports).

**Pipeline (4 stages) — all in scope for this workspace:**
1. **Data collection** — Xsens DOT IMUs, 120 Hz, BLE → thin C++ collector (laptop or
   onboard/shore Raspberry Pi) → raw frames forwarded over TCP/UDP to a server. No
   filtering happens at this stage; keep it thin for real-time. (Live capture
   hardware/networking is lower priority — the dataset below stands in for this
   stage during development.)
2. **Filtering** — adaptive-cutoff Butterworth low-pass filter. **In scope, build
   this.** See section 3 below for the method and section 5 for where the spec lives.
3. **Stroke segmentation & metrics** — zero-crossing detection on the filtered
   Acc_Y, peak detection, per-phase metric extraction. **In scope, build this.**
   See `docs/stroke_detection_spec.md`.
4. **Delivery** — completed-stroke metrics pushed as JSON over WebSocket to browser
   clients, within ~1s of stroke completion. Node.js-based, mobile-compatible web UI.
   Lower priority than 2–3; don't block on it.

**Immediate goal:** implement stages 2 and 3 end-to-end against the attached
dataset (`data/Canoe_xsens_dot.xlsx`) — raw IMU data in, filtered signal, detected
strokes, and per-phase statistics out. Build it causal/streaming-shaped from the
start (see section 3) so it can later run on the live feed from stage 1 without a
rewrite, even though for now it's exercised against a static file.

## 2. Hardware / data facts (confirmed, don't re-derive)

- Sensor: Xsens DOT (Movella), mounted on a crossbar in front of the driver's seat,
  near the boat's center (rider sits slightly aft of center).
- Sample rate: **120 Hz**, confirmed from real data (`SampleTimeFine` deltas = 8333 µs).
- Two sensors per boat (front/rear), each processed independently — do not couple
  their state machines or their filters.
- Channels: `Acc_X, Acc_Y, Acc_Z` (m/s²), `Gyr_X, Gyr_Y, Gyr_Z` (deg/s — confirm
  units against the actual Xsens DOT export before assuming).
- **Acc_Y is the paddling axis.** Sign convention (post-filter): positive =
  acceleration phase (paddle pulling the boat forward), negative = deceleration
  phase (recovery/glide). *Re-verify this against the actual dataset before
  hardcoding* — it depends on sensor mounting orientation.
- Reference/validation file: `Canoe_xsens_dot.xlsx`, sheet `輕艇xsens`, columns
  `Acc_X, Acc_Y, Acc_Z, Gyr_X, Gyr_Y, Gyr_Z` — this is **raw, unfiltered** data.
  It also has a `分期` (phase) column with hand-labeled 0/1 ground truth — use this
  to validate segmentation output, not just to eyeball whether code runs.

## 3. Algorithm decisions

### 3.1 Filtering (now in scope)

Per the project report: a Butterworth low-pass filter (2nd or 4th order — flat
passband, gentle rolloff), but with an **adaptive cutoff** instead of a fixed one,
because a fixed cutoff that works for one paddler/session doesn't generalize (the
report notes 5 Hz worked for one dataset, 10 Hz for another).

**Proposed method (start here, tune against real data — this is not yet validated,
flag results either way):**
1. Maintain a rolling window of raw Acc_Y (2–4s of samples, i.e. 240–480 samples
   @ 120 Hz).
2. Compute a PSD over that window (Welch's method preferred over a raw FFT — more
   robust to noise with short windows). This follows the cutoff-selection
   precedent from PSD-based heart-rate filtering literature (see section 6),
   adapted here for paddling motion instead of cardiac cycles.
3. Find the dominant frequency peak within a plausible paddling-cadence range
   (roughly 0.3–2 Hz, i.e. 20–120 strokes/min) — this is the fundamental stroke
   frequency.
4. Set the cutoff at a multiple of that fundamental (e.g. 4–6x, to preserve enough
   harmonics to keep the stroke's acceleration *shape*, not just its rate), clipped
   to a sane range of **5–20 Hz** (the range repeatedly cited in IMU/human-motion
   filtering literature — see section 6). Don't let the adaptive logic wander
   outside this range; treat an out-of-range result as a bug to investigate, not
   a valid cutoff.
5. Feed the resulting cutoff into a 2nd or 4th order Butterworth low-pass,
   implemented as a causal IIR (cascaded biquads), not a zero-phase/`filtfilt`
   style filter — the real-time system can't look into the future, so the offline
   dataset implementation shouldn't rely on tricks that won't survive the port to
   streaming.
6. Recompute the cutoff periodically as the window slides (e.g. every ~1s) rather
   than every sample — decide how to handle cutoff changes mid-stream without
   introducing filter transients/instability (e.g. only update coefficients at a
   phase boundary, or cross-fade) and document whatever approach is chosen.

If Welch/FFT-based PSD estimation turns out too expensive for the target hardware,
that's worth flagging — don't silently fall back to a fixed cutoff without saying so.

### 3.2 Stroke segmentation

- Zero-crossing on **filtered** Acc_Y, with a **hysteresis band** (start: ±0.15
  m/s²) and a **minimum phase duration** (start: 80 ms / ~10 samples @ 120Hz) to
  reject noise-triggered false crossings. Both must be runtime-configurable, not
  hardcoded — they'll need retuning once real filtered data is available, and
  they interact with whatever the filter stage leaves behind.
- Why not a fancier method (ML classifier, wavelet transform) for segmentation
  itself: threshold/zero-crossing is the accepted first-stage method even in
  papers that add ML on top (e.g. Parab et al. 2025 use thresholding then ML for
  stroke *quality* classification, not segmentation). Don't swap this out for
  something heavier unless real data proves the simple approach insufficient —
  document the failure mode first if you do.

### 3.3 Peak detection

Simple max/min scan within a closed phase's buffer. The signal should already be
clean coming out of the adaptive filter — if real data still shows spurious
sub-peaks after filtering, that's a signal the filter's cutoff selection needs
retuning before reaching for a heavier peak-detection method (see `libPeak`
reference in section 6).

### 3.4 Hard constraints (apply to both filter and segmentation code)

- **Causal only.** Neither the filter nor the segmentation/metrics logic may use
  future samples relative to the point of computation. A phase's metrics are
  computed once that phase has *closed* (confirmed zero-crossing) — that's fine,
  the phase is over, not a look-ahead violation. The filter's rolling-window PSD
  estimate is likewise only ever looking backward.
- **Real-time budget:** stroke JSON must be emitted within ~1s of the stroke's
  deceleration phase closing (filter + segmentation + metrics combined).

## 4. Metric formulas (canonical — see task spec for full tables)

Per closed phase, `dt = 1/120s`, `N` samples in phase:
- `duration_s = N * dt`
- `MAC(x) = mean(|x[i+1] - x[i]|)` (mean absolute change)
- `velocity_integral = trapezoidal_integral(acc_y_phase, dx=dt)`
- `RAD = pos_peak / time_to_peak_seconds` (acceleration phase only)
- Time-to-peak is reported both as a percentage of phase duration and in raw
  seconds (see item 4 / "consolidated" metrics in the task spec) — these are the
  same underlying value, don't compute them independently or they can drift out
  of sync.

Full formula tables, JSON output schema, and pseudocode for the segmentation state
machine live in `docs/stroke_detection_spec.md`. That file's intro currently says
the filter is "owned by a teammate" — that line is stale now that filtering is
in scope here; treat this AGENTS.md as the current source of truth on scope and
flag if the two docs disagree elsewhere.

## 5. Reference materials in this repo

- `docs/stroke_detection_spec.md` — implementation spec for the segmentation/
  metrics module (input struct, state machine pseudocode, formula tables, JSON
  schema, acceptance criteria). Filter section is stale — see section 3 above
  for the current filter approach instead.
- `reference/stroke_metrics_prototype.py` — validated Python reference
  implementation of segmentation + metrics, **currently run on raw (unfiltered)
  data** as a placeholder. Extend this (or its C++ port) to run the adaptive
  filter first, and re-validate the downstream numbers once it does — filtering
  the input will shift peak values, phase durations, etc., so don't assume the
  existing `stroke_metrics.csv` stays valid as ground truth once filtering is added.
- `reference/stroke_metrics.csv` — reference output from the *unfiltered* prototype
  (45 strokes from ~60s of data). Useful for sanity-checking that segmentation
  logic itself still behaves the same shape of way post-filtering, not for exact
  value parity anymore.
- `reference/stroke_segmentation_plot.png`, `reference/stroke_plots/` — visual
  reference for what "correct-looking" segmentation looks like on this dataset
  (pre-filtering) — expect filtered plots to look similar but smoother.
- `data/Canoe_xsens_dot.xlsx` — the raw dataset. This is now the actual input to
  the whole pipeline (filter → segment → metrics), not just a segmentation-stage
  test fixture.

## 6. Literature this project builds on

- Geneau et al. 2024, *Sensors* 24:6085 — zero-crossing rowing stroke segmentation
  (closest methodological precedent for our segmentation approach)
- Dobra et al. 2025, *Sensors* 25:5307 — real-time embedded (ESP32) paddle stroke
  detection architecture (precedent for our "thin collector, heavy server" split)
- Parab et al. 2025, arXiv:2508.01511 — canoe-specific stroke/phase segmentation;
  explicitly notes thresholding reliability issues, which is why we use hysteresis +
  minimum-duration debounce rather than a naive sign check
- Yu et al., "Estimate of the Optimum Cutoff Frequency for the Butterworth
  Low-Pass Digital Filter" — foundational cutoff-selection methodology
- PSD-based cutoff selection for wearable heart-rate filtering (referenced in the
  original project report) — the precedent our adaptive-cutoff method in section
  3.1 is adapted from; that paper used a *fixed* cutoff chosen via PSD, we're
  making it dynamic per rolling window
- General IMU-filtering literature supporting the 5–20 Hz cutoff range used as a
  sanity bound in section 3.1
- Amat et al. 2025, *Algorithms* 18:242 — real-time canoe *training-phase*
  (accel/cruise/sprint) analysis; different granularity than our per-stroke work,
  cite for system-level real-time framing only
- `nerajbobra/embedded_pedometer` (GitHub) — real-time causal event-detection
  algorithm structure, useful analog for the segmentation state machine
- `zmzhang/libPeak` (GitHub) — alternative peak-detection methods if simple
  max/min scanning proves insufficient post-filtering

## 7. Coding conventions for this project

- **Language:** C++ for the real-time server-side pipeline (filter + segmentation +
  delivery); Python is fine for prototyping/validation/offline analysis only —
  don't ship Python in the real-time path.
- **Config over constants:** any tunable threshold (hysteresis, min phase duration,
  filter order, cutoff-selection range/multiplier, etc.) should be a constructor/
  config parameter, not a `#define` or magic number — these all get retuned as
  real data comes in.
- **Bounded memory:** streaming modules (both the filter's PSD window and the
  segmenter's phase buffer) should only hold the data they actually need — never
  buffer a whole session in memory.
- **Two sensors, independent state:** don't share state or add cross-sensor
  coupling between front/rear — for the filter's rolling PSD window as much as
  for the segmenter's phase state — unless a task explicitly asks for
  sensor-fusion work.
- **Validate against ground truth, not just "it compiles":** every filter or
  segmentation change should be checked against `Canoe_xsens_dot.xlsx`'s `分期`
  labels and by re-running the reference prototype, not just run without crashing.
  Since filtering is new, also sanity-check the filter's output directly (e.g.
  plot raw vs. filtered Acc_Y) before trusting anything downstream of it.

## 8. Things NOT to do without asking

- Don't use a single fixed cutoff frequency for the filter — the whole point of
  this stage is that the cutoff adapts per rolling window. If adaptive selection
  isn't working, say so rather than quietly hardcoding a value.
- Don't use non-causal filtering (e.g. `scipy.signal.filtfilt` or any
  forward-backward/zero-phase approach) anywhere in code meant to eventually run
  on the live stream — even in the offline dataset version.
- Don't change the Acc_Y sign convention without first confirming against the
  actual dataset — flag it and ask if something looks inverted, don't silently
  flip it.
- Don't swap the segmentation method (zero-crossing) for an ML-based approach
  without discussing first — see section 3.2 for why.
- Don't invent additional metrics beyond what's in `docs/stroke_detection_spec.md`
  without checking against the original competition requirements first.
- Don't treat `reference/stroke_metrics.csv` as ground truth for exact values once
  filtering is in the loop — it was generated from unfiltered data (see section 5).

## 9. Assumed folder layout (confirm/adjust on first task)

If the actual workspace differs, the first task's agent should confirm real paths
before writing code, not silently assume:

```
project-root/
├── AGENTS.md                              (this file)
├── docs/
│   ├── stroke_detection_spec.md           (segmentation/metrics spec — filter
│   │                                        section is stale, see AGENTS.md §3)
│   └── filter_spec.md                     (create if a dedicated filter spec is
│                                            needed beyond AGENTS.md §3.1)
├── reference/
│   ├── stroke_metrics_prototype.py        (currently unfiltered — extend it)
│   ├── stroke_metrics.csv                 (unfiltered reference output)
│   ├── stroke_segmentation_plot.png
│   └── stroke_plots/
├── data/
│   └── Canoe_xsens_dot.xlsx               (raw data — now the pipeline's actual input)
└── src/                                    (C++ implementation goes here — agent
                                              may need to create this, likely with
                                              filter/ and segmentation/ subfolders)
```

## 10. Status log

Keep a running one-line-per-entry log here as tasks complete, so future agents
(and you) can see what's already built without re-reading full task specs.

- [x] Adaptive-cutoff Butterworth filter (PSD-based, rolling window) — Completed in Python prototype ([src/stroke_analyzer.py](file:///Users/8identity-dev/Workspace/CSA/csa/src/stroke_analyzer.py))
- [x] Stroke segmentation + peak detection + per-phase metrics module (C++) — Completed and verified double-precision parity with Python prototype ([src/StrokeSegmenter.cpp](file:///Users/8identity-dev/Workspace/CSA/csa/src/StrokeSegmenter.cpp))
- [x] Re-validation of segmentation/metrics against filtered (not raw) signal — Completed in Python prototype validation (saved in [reference/stroke_metrics_filtered.csv](file:///Users/8identity-dev/Workspace/CSA/csa/reference/stroke_metrics_filtered.csv))
- [x] Formatted Excel statistics & per-cycle report generator — Completed and integrated into pipeline ([reference/stroke_statistics_report.xlsx](file:///Users/8identity-dev/Workspace/CSA/csa/reference/stroke_statistics_report.xlsx))
