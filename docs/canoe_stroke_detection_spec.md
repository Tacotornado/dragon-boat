# Task: Implement real-time canoe stroke segmentation & metrics module (C++)

## 1. Context

This is one module in a larger dragon boat / canoe racing analytics system. Full pipeline:

1. **Sensor layer** (already exists): Xsens DOT IMUs, 120 Hz, streaming Acc_X/Y/Z + Gyr_X/Y/Z over BLE to a lightweight collector, forwarded raw to a server over TCP/UDP.
2. **Filter stage** (owned by a teammate, already implemented separately): adaptive-cutoff Butterworth low-pass filter (2nd or 4th order), cutoff chosen from a rolling 2–4s PSD/FFT window. Output: a clean Acc_X/Y/Z + Gyr_X/Y/Z stream at the same 120 Hz.
3. **← THIS MODULE:** stroke segmentation (zero-crossing on filtered Acc_Y) + peak detection + per-phase metric extraction, running causally/in real time on the filtered stream.

> **Note:** this spec originally assumed the filter (step 2) was owned separately
> and only consumed here. That's no longer the case — filtering is now also in
> scope for this workspace. See `AGENTS.md` section 3.1 for the adaptive-cutoff
> Butterworth approach to implement upstream of this module. The segmentation
> logic below is unchanged either way; it just now needs to run on a signal your
> own filter code produces, not a teammate's.
4. **Delivery layer** (already exists): completed-stroke metrics get pushed as JSON over WebSocket to browser clients, within ~1 second of stroke completion.

**Your job is step 3 only.** Assume the filtered signal is handed to you sample-by-sample (or in small batches) at 120 Hz; you do not implement the filter itself.

## 2. Input

Per-sample input, 120 Hz, one struct per IMU per sample:

```cpp
struct ImuSample {
    uint64_t timestamp_us;   // monotonic, from sensor
    float acc_x, acc_y, acc_z; // m/s^2, already Butterworth-filtered
    float gyr_x, gyr_y, gyr_z; // deg/s (or rad/s — confirm units with filter teammate), already filtered
};
```

Sensor is mounted on a crossbar near the boat's center, just in front of the driver's seat. **Acc_Y is the paddling axis** (fore-aft), and its sign convention is: positive = acceleration phase (paddle pulling boat forward), negative = deceleration phase (recovery/glide). Confirm this sign convention against real data before hardcoding it — invert if the mounting orientation is different.

A real sample file (for local testing, not the live feed) is attached: `Canoe_xsens_dot.xlsx`, sheet `輕艇xsens`, columns `Acc_X, Acc_Y, Acc_Z, Gyr_X, Gyr_Y, Gyr_Z`, plus a manually-labeled ground-truth column `分期` (0/1, phase label) you can use for validation. Sample rate confirmed at 120 Hz from `SampleTimeFine` deltas (8333 µs).

## 3. Algorithm

### 3.1 Zero-crossing stroke segmentation (causal state machine)

Maintain per-IMU state:

```cpp
enum class Phase { ACCEL, DECEL };

struct SegmenterState {
    Phase current_phase;
    size_t phase_start_idx;       // index (or timestamp) where current phase began
    std::vector<ImuSample> phase_buffer; // samples since phase_start_idx
};
```

Constants (start with these, but expose as tunable config — they may need retuning once real filtered data is available):

```cpp
constexpr float HYSTERESIS = 0.15f;      // m/s^2 dead-band around zero, avoids noise-triggered crossings
constexpr float MIN_PHASE_DURATION_S = 0.08f; // reject phases shorter than this as noise (~10 samples @ 120Hz)
```

Per-sample update logic:

```
on new sample s:
    append s to phase_buffer
    if current_phase == ACCEL and s.acc_y < -HYSTERESIS:
        if (s.timestamp - phase_start_time) >= MIN_PHASE_DURATION_S:
            emit_completed_phase(ACCEL, phase_buffer)   // triggers metric calc, see 3.3
            if previous emitted phase was DECEL (i.e. we just completed a full stroke):
                emit_completed_stroke(...)               // combine last ACCEL+DECEL pair
            phase_start_time = s.timestamp
            phase_buffer = [s]
            current_phase = DECEL
    elif current_phase == DECEL and s.acc_y > HYSTERESIS:
        if (s.timestamp - phase_start_time) >= MIN_PHASE_DURATION_S:
            emit_completed_phase(DECEL, phase_buffer)
            phase_start_time = s.timestamp
            phase_buffer = [s]
            current_phase = ACCEL
    # else: still inside current phase, keep buffering
```

Important: this is a **causal, streaming** algorithm — never look ahead. The only latency is the buffer for the phase currently in progress, which closes naturally when the next zero-crossing is confirmed. A stroke's metrics become available immediately after its deceleration phase closes, satisfying the "~1s after stroke completion" delivery requirement.

A short false-positive phase (e.g., noise causing state to flip and flip back within a few samples) should be rejected by `MIN_PHASE_DURATION_S` — verify this against real data, since even filtered signals can wobble near zero on a rough water surface.

### 3.2 Peak detection

Within each completed phase's buffer, find:
- **Acceleration phase**: `pos_peak = max(acc_y)` over the buffer, and its index/timestamp.
- **Deceleration phase**: `neg_peak = min(acc_y)` over the buffer, and its index/timestamp.

This is a simple linear-scan max/min — no smoothing filter needed here since the signal is already Butterworth-filtered upstream. If real data shows spurious sub-peaks even after filtering, consider a minimum-prominence check (reject a peak if it's within X% of the phase's other local extrema) before falling back to more complex methods (wavelet/derivative-based, see `libPeak` reference below).

### 3.3 Per-phase metrics (compute at the moment each phase closes)

Let `dt = 1/120` s. For a phase spanning samples `[s, e]` (start index `s`, end index `e`, buffer of `N = e - s` samples):

```
duration_s          = N * dt
mean_absolute_change(x) = mean(|x[i+1] - x[i]|) for i in [s, e-1]      // MAC
velocity_integral   = trapezoidal_integral(acc_y[s:e], dx=dt)          // "change in velocity"
```

**Acceleration phase, on close:**
| Metric | Formula |
|---|---|
| Acceleration phase duration | `duration_s` |
| ACC-Y positive peak | `pos_peak` (value) |
| Time to ACC-Y positive peak (%) | `100 * (t_peak - t_start) / duration_s` |
| RAD (rate of acceleration development) | `pos_peak / (t_peak - t_start)` — guard against divide-by-zero if peak is at index 0 |
| Velocity change (integral) | `velocity_integral` over this phase |
| MAC of ACC-X | `mean_absolute_change(acc_x)` over this phase |
| MAC of Gyr-X, Gyr-Y, Gyr-Z | same, per axis |

**Deceleration phase, on close:**
| Metric | Formula |
|---|---|
| Deceleration phase duration | `duration_s` |
| ACC-Y negative peak | `neg_peak` (value) |
| Time to ACC-Y negative peak (%) | `100 * (t_peak - t_start) / duration_s` |
| Velocity change (integral) | `velocity_integral` over this phase |
| MAC of ACC-X, Gyr-X, Gyr-Y, Gyr-Z | same as above, over this phase |

**Whole stroke (accel phase immediately followed by decel phase — emit once decel closes):**
| Metric | Formula |
|---|---|
| Gyr-X positive peak | `max(gyr_x)` over accel_start → decel_end |
| Gyr-X negative peak | `min(gyr_x)` over accel_start → decel_end |
| MAC of ACC-X, Gyr-X, Gyr-Y, Gyr-Z | over the full accel_start → decel_end span |

**Consolidated raw-time summary (item 4 in the original spec — equivalent quantities, non-percentage):**
| Metric | Formula |
|---|---|
| TotalTime | `accel_duration_s + decel_duration_s` |
| TimeToPosPeak | `t_pos_peak - t_accel_start` (seconds, not %) |
| TimeToNegPeak | `t_neg_peak - t_decel_start` (seconds, not %; "mid zero" = the accel→decel crossing) |
| Area_Positive | same as the acceleration phase's velocity_integral |
| Area_Negative | same as the deceleration phase's velocity_integral |

A full worked Python reference implementation (validated against the sample `.xlsx` data, 45 strokes detected from ~60s) is attached as `stroke_metrics_prototype.py` — port its logic directly rather than re-deriving it; it already handles edge cases like the first/last incomplete phase in a session.

## 4. Output

Emit one JSON object per completed stroke (accel phase + its following decel phase):

```json
{
  "stroke_id": 7,
  "sensor_id": "front",
  "accel": {
    "duration_s": 0.425,
    "acc_y_pos_peak": 5.781,
    "time_to_peak_pct": 37.25,
    "rad": 36.51,
    "velocity_integral": 1.423,
    "mac_acc_x": 0.110, "mac_gyr_x": 0.273, "mac_gyr_y": 1.044, "mac_gyr_z": 0.068
  },
  "decel": {
    "duration_s": 1.117,
    "acc_y_neg_peak": -2.156,
    "time_to_peak_pct": 81.34,
    "velocity_integral": -1.561,
    "mac_acc_x": 0.123, "mac_gyr_x": 0.169, "mac_gyr_y": 1.415, "mac_gyr_z": 0.035
  },
  "stroke": {
    "gyr_x_pos_peak": 7.564, "gyr_x_neg_peak": -4.063,
    "mac_acc_x": 0.119, "mac_gyr_x": 0.197, "mac_gyr_y": 1.320, "mac_gyr_z": 0.044,
    "total_time_s": 1.542, "time_to_pos_peak_s": 0.158, "time_to_neg_peak_s": 0.908,
    "area_positive": 1.423, "area_negative": -1.561
  }
}
```

This gets handed to the existing WebSocket delivery layer — do not implement transport here, just produce this struct/JSON.

## 5. Non-functional requirements

- **Causal only.** No processing may use future samples relative to the point of computation, except within a phase's own buffer once that phase has closed (this is legitimate — the phase is over, not the future).
- **Per-IMU independence.** Two sensors (front/rear) run this pipeline independently; don't couple their state.
- **Bounded memory.** Only the currently-open phase's samples need to be buffered (typically well under 1s of data = ~120 samples); nothing session-long should be kept in memory by this module.
- **Config, not constants.** `HYSTERESIS` and `MIN_PHASE_DURATION_S` should be runtime-configurable (struct/constructor params), since they will need retuning once real filtered (not raw) data is available.

## 6. Validation plan

1. Unit tests with synthetic data: a clean sine wave on Acc_Y should produce exactly N/2 accel + N/2 decel phases of predictable duration; inject noise near zero-crossings and confirm the hysteresis + min-duration logic suppresses false phases.
2. Run against the attached `Canoe_xsens_dot.xlsx` (`輕艇xsens` sheet) and compare detected phase boundaries against the manually-labeled `分期` column — report agreement (e.g. overlap/IoU per phase, or simple boundary-timestamp delta).
3. Cross-check a handful of strokes' computed metrics against the attached `stroke_metrics.csv` (reference output from the Python prototype) to confirm the C++ port produces matching numbers within floating-point tolerance.
4. Visual spot-check: reproduce the per-stroke plots (see attached `stroke_plots/` folder) from the C++ output to confirm peaks and zero-crossings line up with what a human would call a "clean" stroke.

## 7. Reference materials (attach alongside this spec)

- `stroke_metrics_prototype.py` — working Python reference implementation + validated formulas
- `stroke_metrics.csv` — reference output (45 strokes) to check numeric parity against
- `stroke_segmentation_plot.png` — overview plot (zoomed + full session) showing the expected segmentation visually
- `stroke_plots/` — one annotated plot per detected stroke (full cycle, peaks + zero-crossings marked)
- `Canoe_xsens_dot.xlsx` — real sample data with ground-truth phase labels for validation

Relevant papers for methodology justification (not required reading to implement, useful for the competition report):
- Geneau et al. 2024, *Sensors* 24:6085 — zero-crossing-based rowing stroke phase segmentation (same technique family)
- Dobra et al. 2025, *Sensors* 25:5307 — real-time embedded (ESP32) paddle stroke detection architecture
- Parab et al. 2025 (arXiv:2508.01511) — canoe-specific stroke/phase segmentation pipeline, notes threshold-based segmentation reliability issues worth guarding against
- Yu et al. — Butterworth cutoff selection methodology (already in use by the filter teammate)

## 8. Acceptance criteria

- [ ] Compiles as a standalone C++ module with a clean interface (`process_sample(ImuSample) -> optional<StrokeResult>`)
- [ ] Produces stroke count and per-stroke metrics matching the Python prototype's CSV within reasonable floating-point tolerance, on the same input file
- [ ] Correctly rejects short noise-triggered phases (verify with the synthetic sine + noise unit test)
- [ ] No look-ahead / non-causal logic anywhere in the hot path
- [ ] Runtime-configurable hysteresis and minimum phase duration
- [ ] Emits stroke JSON within the existing "~1s after stroke completion" delivery budget (verify with a timestamp check between decel-phase-close and JSON-emit)
