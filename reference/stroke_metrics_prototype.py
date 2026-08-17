"""
Canoe stroke segmentation + metric extraction prototype.

Pipeline stage this script covers (per the project report):
  Input  : filtered Acc/Gyr streams (here: raw front-sensor data, since the
           Butterworth-filtered signal from your teammate isn't in this file yet -
           swap `acc_y_raw` etc. for the filtered arrays once you have them).
  Output : one row per completed stroke with all Acceleration / Deceleration /
           Whole-stroke metrics requested in the spec.

This is written as an *offline* prototype (reads the whole file, uses numpy/pandas)
so you can validate the math fast. Section 2 shows how the same logic maps onto
a causal, sample-by-sample state machine for the real-time C++ port.
"""

import numpy as np
import pandas as pd

FS = 120.0          # Hz, Xsens DOT sample rate
DT = 1.0 / FS

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
df = pd.read_excel(
    "/mnt/user-data/uploads/Canoe_xsens_dot.xlsx",
    sheet_name="輕艇xsens",
    usecols="A:H",
)
df = df.dropna(subset=["Acc_Y"]).reset_index(drop=True)

acc_x = df["Acc_X"].to_numpy()
acc_y = df["Acc_Y"].to_numpy()
gyr_x = df["Gyr_X"].to_numpy()
gyr_y = df["Gyr_Y"].to_numpy()
gyr_z = df["Gyr_Z"].to_numpy()
n = len(acc_y)

# ---------------------------------------------------------------------------
# 2. Zero-crossing segmentation on ACC_Y, with hysteresis + debounce
#    (this is the part that maps directly onto a real-time state machine)
# ---------------------------------------------------------------------------
HYST = 0.15          # m/s^2 dead-band around zero -> ignores small noise wobble
MIN_PHASE_SAMPLES = int(0.08 * FS)   # ignore phases shorter than 80 ms (noise)

state = "accel" if acc_y[0] >= 0 else "decel"
phase_start = 0
phases = []   # list of dicts: {kind, start, end}

for i in range(1, n):
    if state == "accel" and acc_y[i] < -HYST:
        if i - phase_start >= MIN_PHASE_SAMPLES:
            phases.append({"kind": "accel", "start": phase_start, "end": i})
            phase_start = i
            state = "decel"
    elif state == "decel" and acc_y[i] > HYST:
        if i - phase_start >= MIN_PHASE_SAMPLES:
            phases.append({"kind": "decel", "start": phase_start, "end": i})
            phase_start = i
            state = "accel"
# close the last open phase
phases.append({"kind": state, "start": phase_start, "end": n - 1})

# ---------------------------------------------------------------------------
# 3. Per-phase metrics
# ---------------------------------------------------------------------------
def mac(x):
    """Mean absolute change: mean(|x[i+1]-x[i]|)"""
    return np.mean(np.abs(np.diff(x))) if len(x) > 1 else 0.0

rows = []
i = 0
while i < len(phases) - 1:
    p = phases[i]
    if p["kind"] != "accel":
        i += 1
        continue
    # need an accel phase immediately followed by a decel phase to form one stroke
    if i + 1 >= len(phases) or phases[i + 1]["kind"] != "decel":
        i += 1
        continue

    dphase = phases[i + 1]
    a_s, a_e = p["start"], p["end"]
    d_s, d_e = dphase["start"], dphase["end"]

    seg_ay_a = acc_y[a_s:a_e]
    seg_ay_d = acc_y[d_s:d_e]

    accel_dur = (a_e - a_s) * DT
    decel_dur = (d_e - d_s) * DT

    pos_peak_idx = a_s + int(np.argmax(seg_ay_a))
    neg_peak_idx = d_s + int(np.argmin(seg_ay_d))

    t_to_pos_peak = (pos_peak_idx - a_s) * DT
    t_to_neg_peak = (neg_peak_idx - d_s) * DT

    pos_peak = acc_y[pos_peak_idx]
    neg_peak = acc_y[neg_peak_idx]

    row = {
        "stroke_id": len(rows) + 1,
        # --- Acceleration phase ---
        "accel_duration_s": accel_dur,
        "accY_pos_peak": pos_peak,
        "t_to_pos_peak_pct": 100 * t_to_pos_peak / accel_dur if accel_dur > 0 else np.nan,
        "RAD": pos_peak / t_to_pos_peak if t_to_pos_peak > 0 else np.nan,
        "accel_velocity_integral": np.trapezoid(seg_ay_a, dx=DT),
        "accel_MAC_accX": mac(acc_x[a_s:a_e]),
        "accel_MAC_gyrX": mac(gyr_x[a_s:a_e]),
        "accel_MAC_gyrY": mac(gyr_y[a_s:a_e]),
        "accel_MAC_gyrZ": mac(gyr_z[a_s:a_e]),
        # --- Deceleration phase ---
        "decel_duration_s": decel_dur,
        "accY_neg_peak": neg_peak,
        "t_to_neg_peak_pct": 100 * t_to_neg_peak / decel_dur if decel_dur > 0 else np.nan,
        "decel_velocity_integral": np.trapezoid(seg_ay_d, dx=DT),
        "decel_MAC_accX": mac(acc_x[d_s:d_e]),
        "decel_MAC_gyrX": mac(gyr_x[d_s:d_e]),
        "decel_MAC_gyrY": mac(gyr_y[d_s:d_e]),
        "decel_MAC_gyrZ": mac(gyr_z[d_s:d_e]),
        # --- Whole stroke (accel + decel) ---
        "stroke_gyrX_pos_peak": np.max(gyr_x[a_s:d_e]),
        "stroke_gyrX_neg_peak": np.min(gyr_x[a_s:d_e]),
        "stroke_MAC_accX": mac(acc_x[a_s:d_e]),
        "stroke_MAC_gyrX": mac(gyr_x[a_s:d_e]),
        "stroke_MAC_gyrY": mac(gyr_y[a_s:d_e]),
        "stroke_MAC_gyrZ": mac(gyr_z[a_s:d_e]),
        # --- Item 4: consolidated raw-time summary ---
        "TotalTime_s": accel_dur + decel_dur,
        "TimeToPosPeak_s": t_to_pos_peak,
        "TimeToNegPeak_s": t_to_neg_peak,
        "Area_Positive": np.trapezoid(seg_ay_a, dx=DT),
        "Area_Negative": np.trapezoid(seg_ay_d, dx=DT),
    }
    rows.append(row)
    i += 2

result = pd.DataFrame(rows)
result.to_csv("/mnt/user-data/outputs/stroke_metrics.csv", index=False)
print(f"Detected {len(result)} strokes from {n} samples ({n/FS:.1f} s of data)")
print(result.head(10).to_string())
