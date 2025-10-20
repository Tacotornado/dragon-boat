# stroke_cycle_dtt_latest.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, butter, filtfilt
import os

# --- Setup ---
RESULT_DIR = "results_full_strokes"
STROKE_PLOT_DIR = os.path.join(RESULT_DIR, "stroke_plots")
os.makedirs(STROKE_PLOT_DIR, exist_ok=True)

# --- Load Data ---
df = pd.read_excel("data/Canoe xsens dot.xlsx", sheet_name="輕艇xsens")
acc_y = df['Acc_Y'].values

# --- Filtering ---
def butter_lowpass(cutoff, fs, order=4):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def lowpass_filter(data, cutoff=10, fs=120, order=4):
    b, a = butter_lowpass(cutoff, fs, order=order)
    return filtfilt(b, a, data)

fs = 120
acc_y_smooth = lowpass_filter(acc_y, cutoff=5, fs=fs, order=4)

# --- Zero crossing detection ---
def filter_close_zero_crossings(zcrossings, acc, window=7, min_separation=10):
    filtered = []
    zcrossings = sorted(zcrossings)
    i = 0
    while i < len(zcrossings):
        current = zcrossings[i]
        cluster = [current]
        j = i + 1
        while j < len(zcrossings) and abs(zcrossings[j] - current) < min_separation:
            cluster.append(zcrossings[j])
            j += 1
        if len(cluster) == 1:
            filtered.append(current)
        else:
            amps = []
            for z in cluster:
                zero_idx = int(round(z))
                lo = max(0, zero_idx - window)
                hi = min(len(acc), zero_idx + window + 1)
                local_segment = acc[lo:hi]
                local_amp = np.max(local_segment) - np.min(local_segment)
                amps.append(local_amp)
            max_idx = np.argmax(amps)
            filtered.append(cluster[max_idx])
        i += len(cluster)
    return np.array(filtered)

def find_all_zero_crossings(acc, min_amplitude=0.5, window=7, eps=0.05):
    zero_crossings = []
    n = len(acc)
    i = 0
    while i < n - 1:
        if (acc[i] > eps and acc[i+1] < -eps) or (acc[i] < -eps and acc[i+1] > eps):
            denom = acc[i+1] - acc[i]
            z = i - acc[i] / denom
            zero_idx = int(round(z))
            lo = max(0, zero_idx - window)
            hi = min(n, zero_idx + window + 1)
            local_segment = acc[lo:hi]
            local_amp = np.max(local_segment) - np.min(local_segment)
            slope = abs(acc[i+1] - acc[i])
            if local_amp >= min_amplitude and slope > 0.15:
                zero_crossings.append(z)
            i += 1
        elif abs(acc[i]) <= eps:
            start = i
            while i < n and abs(acc[i]) <= eps:
                i += 1
            end = i - 1
            mid = (start + end) / 2.0
            zero_idx = int(round(mid))
            lo = max(0, zero_idx - window)
            hi = min(n, zero_idx + window + 1)
            local_segment = acc[lo:hi]
            local_amp = np.max(local_segment) - np.min(local_segment)
            if local_amp >= min_amplitude:
                zero_crossings.append(mid)
        else:
            i += 1
    return filter_close_zero_crossings(zero_crossings, acc, window)

# --- Peaks ---
pos_peaks, _ = find_peaks(acc_y_smooth, distance=5, prominence=0.8)
neg_peaks, _ = find_peaks(-acc_y_smooth, distance=5, prominence=0.8)

# --- Zero crossings ---
all_zeros = find_all_zero_crossings(acc_y_smooth, min_amplitude=0.65, window=7, eps=0.1)
all_zeros_filtered = np.array([z for z in all_zeros if (
    np.max(acc_y_smooth[int(z)-7:int(z)+8]) - np.min(acc_y_smooth[int(z)-7:int(z)+8])
) > 0.65])

# --- Full stroke cycles (0 → +peak → 0 → -peak → 0) ---
stroke_cycles = []
for i in range(len(all_zeros_filtered) - 3):
    z1 = int(all_zeros_filtered[i])
    z2 = int(all_zeros_filtered[i+1])
    z3 = int(all_zeros_filtered[i+2])
    z4 = int(all_zeros_filtered[i+3])

    pos_in_seg = pos_peaks[(pos_peaks > z1) & (pos_peaks < z2)]
    neg_in_seg = neg_peaks[(neg_peaks > z2) & (neg_peaks < z3)]

    if len(pos_in_seg) == 0 or len(neg_in_seg) == 0:
        continue

    pos_idx = pos_in_seg[np.argmax(acc_y_smooth[pos_in_seg])]
    neg_idx = neg_in_seg[np.argmin(acc_y_smooth[neg_in_seg])]

    t_total = (z3 - z1) / fs
    t_pos = (pos_idx - z1) / fs
    t_neg = (neg_idx - z2) / fs

    pos_area = np.trapz(np.clip(acc_y_smooth[z1:z2+1], 0, None), dx=1/fs)
    neg_area = np.trapz(np.clip(acc_y_smooth[z2:z3+1], None, 0), dx=1/fs)

    stroke_cycles.append({
        "StrokeID": len(stroke_cycles) + 1,
        "StartZero": z1,
        "PosPeak": pos_idx,
        "MidZero": z2,
        "NegPeak": neg_idx,
        "EndZero": z3,
        "TotalTime": t_total,
        "TimeToPosPeak": t_pos,
        "TimeToNegPeak": t_neg,
        "Area_Positive": pos_area,
        "Area_Negative": neg_area
    })

# --- Save stroke cycles ---
stroke_df = pd.DataFrame(stroke_cycles)
stroke_df.to_csv(os.path.join(RESULT_DIR, "stroke_full_cycles.csv"), index=False)
print(f"Saved {len(stroke_cycles)} full stroke cycles")

# --- Plot overview ---
plt.figure(figsize=(16, 6))
plt.plot(acc_y_smooth, label='Smoothed Acc_Y', linewidth=2)

pos_cycle_peaks = [sc["PosPeak"] for sc in stroke_cycles]
neg_cycle_peaks = [sc["NegPeak"] for sc in stroke_cycles]
zeros_cycle = []
for sc in stroke_cycles:
    zeros_cycle.extend([sc["StartZero"], sc["MidZero"], sc["EndZero"]])

plt.plot(pos_cycle_peaks, acc_y_smooth[pos_cycle_peaks], 'go', label='Positive Peaks')
plt.plot(neg_cycle_peaks, acc_y_smooth[neg_cycle_peaks], 'ro', label='Negative Peaks')
plt.plot(zeros_cycle, np.zeros(len(zeros_cycle)), 'mo', label='Zero Crossings')

for sc in stroke_cycles:
    x_fill = np.arange(sc["StartZero"], sc["EndZero"]+1)
    y_fill = acc_y_smooth[x_fill]
    plt.fill_between(x_fill, 0, y_fill, color='gray', alpha=0.3)

plt.title("Full Stroke Cycles: 0 → +Peak → 0 → -Peak → 0")
plt.xlabel("Sample Index")
plt.ylabel("Acc_Y")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(RESULT_DIR, "strokes_full_overview.png"), dpi=300)
plt.close()

# --- Individual stroke plots ---
for i, sc in enumerate(stroke_cycles):
    seg_x = np.arange(sc["StartZero"], sc["EndZero"]+1)
    seg_y = acc_y_smooth[seg_x]

    plt.figure(figsize=(10, 4))
    plt.plot(seg_x, seg_y, label="Smoothed Acc_Y", linewidth=2)
    plt.axhline(0, color='k', linestyle='--')
    plt.plot(sc["PosPeak"], acc_y_smooth[sc["PosPeak"]], 'go', label="Positive Peak")
    plt.plot(sc["NegPeak"], acc_y_smooth[sc["NegPeak"]], 'ro', label="Negative Peak")
    plt.plot([sc["StartZero"], sc["MidZero"], sc["EndZero"]], [0,0,0], 'mo', label="Zero Crossings")

    plt.title(f"Full Stroke {i+1}: 0→+Peak→0→-Peak→0")
    plt.xlabel("Sample Index")
    plt.ylabel("Acc_Y")
    plt.legend()
    plt.grid(True)
    seg_path = os.path.join(STROKE_PLOT_DIR, f"stroke_full_{i+1}.png")
    plt.savefig(seg_path, dpi=300)
    plt.close()

print(f"Saved {len(stroke_cycles)} individual stroke plots in '{STROKE_PLOT_DIR}/'")

# --- Phase Metrics (Recalculated per full stroke) ---
acc_x = df['Acc_X'].values
gyr_x = df['Gyr_X'].values
gyr_y = df['Gyr_Y'].values
gyr_z = df['Gyr_Z'].values

acc_x_smooth = lowpass_filter(acc_x, cutoff=5, fs=120, order=4)
gyr_x_smooth = lowpass_filter(gyr_x, cutoff=5, fs=120, order=4)
gyr_y_smooth = lowpass_filter(gyr_y, cutoff=5, fs=120, order=4)
gyr_z_smooth = lowpass_filter(gyr_z, cutoff=5, fs=120, order=4)

phase_metrics = []
for sc in stroke_cycles:
    start_idx = sc["StartZero"]
    mid_idx   = sc["MidZero"]
    end_idx   = sc["EndZero"]
    pos_idx   = sc["PosPeak"]
    neg_idx   = sc["NegPeak"]

    # --- Acceleration phase (StartZero → PosPeak)
    acc_phase = acc_y_smooth[start_idx:pos_idx+1]
    dur_acc = (pos_idx - start_idx) / fs
    acc_peak = np.max(acc_phase)
    time_to_peak = dur_acc
    perc_time_to_peak = 100.0
    RAD = acc_peak / time_to_peak if time_to_peak > 0 else np.nan
    vel_change_acc = np.trapz(acc_phase, dx=1/fs)

    mean_abs_ax_acc  = np.mean(np.abs(np.diff(acc_x_smooth[start_idx:pos_idx+1])))
    mean_abs_gx_acc  = np.mean(np.abs(np.diff(gyr_x_smooth[start_idx:pos_idx+1])))
    mean_abs_gy_acc  = np.mean(np.abs(np.diff(gyr_y_smooth[start_idx:pos_idx+1])))
    mean_abs_gz_acc  = np.mean(np.abs(np.diff(gyr_z_smooth[start_idx:pos_idx+1])))

    # --- Deceleration phase (PosPeak → EndZero)
    dec_phase = acc_y_smooth[pos_idx:end_idx+1]
    dur_dec = (end_idx - pos_idx) / fs
    dec_peak = np.min(dec_phase)
    time_to_neg_peak = (neg_idx - pos_idx) / fs
    perc_time_to_neg = 100 * time_to_neg_peak / dur_dec if dur_dec > 0 else np.nan
    vel_change_dec = np.trapz(dec_phase, dx=1/fs)

    mean_abs_ax_dec  = np.mean(np.abs(np.diff(acc_x_smooth[pos_idx:end_idx+1])))
    mean_abs_gx_dec  = np.mean(np.abs(np.diff(gyr_x_smooth[pos_idx:end_idx+1])))
    mean_abs_gy_dec  = np.mean(np.abs(np.diff(gyr_y_smooth[pos_idx:end_idx+1])))
    mean_abs_gz_dec  = np.mean(np.abs(np.diff(gyr_z_smooth[pos_idx:end_idx+1])))

    # --- Stroke phase (whole cycle StartZero → EndZero)
    gyrx_phase = gyr_x_smooth[start_idx:end_idx+1]
    gyrx_pos_peak = np.max(gyrx_phase)
    gyrx_neg_peak = np.min(gyrx_phase)

    mean_abs_ax_stroke = np.mean(np.abs(np.diff(acc_x_smooth[start_idx:end_idx+1])))
    mean_abs_gx_stroke = np.mean(np.abs(np.diff(gyr_x_smooth[start_idx:end_idx+1])))
    mean_abs_gy_stroke = np.mean(np.abs(np.diff(gyr_y_smooth[start_idx:end_idx+1])))
    mean_abs_gz_stroke = np.mean(np.abs(np.diff(gyr_z_smooth[start_idx:end_idx+1])))

    phase_metrics.append({
            "StrokeID": sc["StrokeID"],
            # Acceleration
            "AccPhase_Duration": dur_acc,
            "AccPhase_PosPeak": acc_peak,
            "AccPhase_TimeToPosPeak_%": perc_time_to_peak,
            "AccPhase_RAD": RAD,
            "AccPhase_VelChange": vel_change_acc,
            "AccPhase_MeanAbs_AccX": mean_abs_ax_acc,
            "AccPhase_MeanAbs_GyrX": mean_abs_gx_acc,
            "AccPhase_MeanAbs_GyrY": mean_abs_gy_acc,
            "AccPhase_MeanAbs_GyrZ": mean_abs_gz_acc,
            # Deceleration
            "DecPhase_Duration": dur_dec,
            "DecPhase_NegPeak": dec_peak,
            "DecPhase_TimeToNegPeak_%": perc_time_to_neg,
            "DecPhase_VelChange": vel_change_dec,
            "DecPhase_MeanAbs_AccX": mean_abs_ax_dec,
            "DecPhase_MeanAbs_GyrX": mean_abs_gx_dec,
            "DecPhase_MeanAbs_GyrY": mean_abs_gy_dec,
            "DecPhase_MeanAbs_GyrZ": mean_abs_gz_dec,
            # Stroke
            "Stroke_GyrX_PosPeak": gyrx_pos_peak,
            "Stroke_GyrX_NegPeak": gyrx_neg_peak,
            "Stroke_MeanAbs_AccX": mean_abs_ax_stroke,
            "Stroke_MeanAbs_GyrX": mean_abs_gx_stroke,
            "Stroke_MeanAbs_GyrY": mean_abs_gy_stroke,
            "Stroke_MeanAbs_GyrZ": mean_abs_gz_stroke,
        })
        
# Convert to DataFrame
phase_df = pd.DataFrame(phase_metrics)

# Split into 3 DataFrames
acc_phase_df = phase_df[[
    "StrokeID",
    "AccPhase_Duration",
    "AccPhase_PosPeak",
    "AccPhase_TimeToPosPeak_%",
    "AccPhase_RAD",
    "AccPhase_VelChange",
    "AccPhase_MeanAbs_AccX",
    "AccPhase_MeanAbs_GyrX",
    "AccPhase_MeanAbs_GyrY",
    "AccPhase_MeanAbs_GyrZ"
]]

dec_phase_df = phase_df[[
    "StrokeID",
    "DecPhase_Duration",
    "DecPhase_NegPeak",
    "DecPhase_TimeToNegPeak_%",
    "DecPhase_VelChange",
    "DecPhase_MeanAbs_AccX",
    "DecPhase_MeanAbs_GyrX",
    "DecPhase_MeanAbs_GyrY",
    "DecPhase_MeanAbs_GyrZ"
]]

stroke_phase_df = phase_df[[
    "StrokeID",
    "Stroke_GyrX_PosPeak",
    "Stroke_GyrX_NegPeak",
    "Stroke_MeanAbs_AccX",
    "Stroke_MeanAbs_GyrX",
    "Stroke_MeanAbs_GyrY",
    "Stroke_MeanAbs_GyrZ"
]]

# Save to one Excel with 3 sheets
phase_excel = os.path.join(RESULT_DIR, "stroke_phase_metrics.xlsx")
with pd.ExcelWriter(phase_excel, engine="openpyxl") as writer:
    acc_phase_df.to_excel(writer, sheet_name="AccelerationPhase", index=False)
    dec_phase_df.to_excel(writer, sheet_name="DecelerationPhase", index=False)
    stroke_phase_df.to_excel(writer, sheet_name="StrokePhase", index=False)

print(f"Phase metrics saved to '{phase_excel}' with 3 sheets")