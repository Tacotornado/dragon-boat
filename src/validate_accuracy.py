import os
import pandas as pd
import numpy as np

def main():
    workspace_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    excel_path = os.path.join(workspace_root, "data/Canoe xsens dot.xlsx")
    filtered_path = os.path.join(workspace_root, "reference/filtered_imu_data.csv")
    metrics_path = os.path.join(workspace_root, "reference/stroke_metrics_cpp.csv")

    if not os.path.exists(excel_path):
        print(f"Error: Excel file not found at {excel_path}")
        return
    if not os.path.exists(filtered_path) or not os.path.exists(metrics_path):
        print("Error: Filtered data or stroke metrics not found in reference/. Please run C++ engine first.")
        return

    print("Loading datasets for validation...")
    df = pd.read_excel(excel_path, sheet_name="輕艇xsens").dropna(subset=["Acc_Y"]).reset_index(drop=True)
    ground_truth = df["分期"].to_numpy()

    filt_df = pd.read_csv(filtered_path)
    acc_y_filt = filt_df["Acc_Y"].to_numpy()
    n = len(acc_y_filt)

    FS = 120.0
    HYST = 0.15
    MIN_PHASE_SAMPLES = 10

    # Reconstruct detected labels for sample-by-sample matching
    our_labels = np.zeros(n)
    state = 1 if acc_y_filt[0] >= 0 else 0
    phase_start = 0
    for i in range(1, n):
        if state == 1 and acc_y_filt[i] < -HYST:
            if i - phase_start >= MIN_PHASE_SAMPLES:
                our_labels[phase_start:i] = 1
                phase_start = i
                state = 0
        elif state == 0 and acc_y_filt[i] > HYST:
            if i - phase_start >= MIN_PHASE_SAMPLES:
                our_labels[phase_start:i] = 0
                phase_start = i
                state = 1
    our_labels[phase_start:n] = state

    # 1. Classification Error
    raw_error = (1.0 - np.mean(our_labels == ground_truth)) * 100
    raw_accuracy = 100.0 - raw_error

    # 2. Shifted Error (Compensating for the 9-sample phase lag)
    shifted_labels = np.zeros(n)
    shift = 9
    shifted_labels[:-shift] = our_labels[shift:]
    shifted_labels[-shift:] = our_labels[-1]

    shifted_error = (1.0 - np.mean(shifted_labels == ground_truth)) * 100
    shifted_accuracy = 100.0 - shifted_error

    # 3. Match strokes and calculate duration percent error
    # Ground truth phases
    gt_changes = np.where(ground_truth[1:] != ground_truth[:-1])[0] + 1
    gt_phases = []
    start = 0
    for change in gt_changes:
        gt_phases.append({"kind": "accel" if ground_truth[start] == 1 else "decel", "start": start, "end": change})
        start = change
    gt_phases.append({"kind": "accel" if ground_truth[start] == 1 else "decel", "start": start, "end": n})

    # Filter out complete strokes in Ground Truth
    gt_strokes = []
    i = 0
    while i < len(gt_phases) - 1:
        if gt_phases[i]["kind"] == "accel" and gt_phases[i+1]["kind"] == "decel":
            gt_strokes.append({
                "accel_len": gt_phases[i]["end"] - gt_phases[i]["start"],
                "decel_len": gt_phases[i+1]["end"] - gt_phases[i+1]["start"],
                "total_len": gt_phases[i+1]["end"] - gt_phases[i]["start"]
            })
            i += 2
        else:
            i += 1

    # Loaded C++ metrics
    cpp_metrics = pd.read_csv(metrics_path)
    detected_strokes = []
    for idx, row in cpp_metrics.iterrows():
        detected_strokes.append({
            "accel_len": int(round(row["accel_duration_s"] * FS)),
            "decel_len": int(round(row["decel_duration_s"] * FS)),
            "total_len": int(round(row["TotalTime_s"] * FS))
        })

    num_strokes_to_compare = min(len(gt_strokes), len(detected_strokes))
    
    accel_errs = []
    decel_errs = []
    total_errs = []

    for idx in range(num_strokes_to_compare):
        gt = gt_strokes[idx]
        det = detected_strokes[idx]
        
        acc_err = abs(det["accel_len"] - gt["accel_len"]) / gt["accel_len"] * 100
        dec_err = abs(det["decel_len"] - gt["decel_len"]) / gt["decel_len"] * 100
        tot_err = abs(det["total_len"] - gt["total_len"]) / gt["total_len"] * 100
        
        accel_errs.append(acc_err)
        decel_errs.append(dec_err)
        total_errs.append(tot_err)

    mean_acc_err = np.mean(accel_errs)
    mean_dec_err = np.mean(decel_errs)
    mean_tot_err = np.mean(total_errs)

    print("\n==================================================")
    print("           VALIDATION SUMMARY REPORT             ")
    print("==================================================")
    print(f"Total Reference Ground Truth Strokes: {len(gt_strokes)}")
    print(f"Total Pipeline Detected Strokes     : {len(detected_strokes)}")
    print("--------------------------------------------------")
    print(f"Raw Sample Accuracy                 : {raw_accuracy:.2f}% (Error: {raw_error:.2f}%)")
    print(f"Lag-Compensated Sample Accuracy     : {shifted_accuracy:.2f}% (Error: {shifted_error:.2f}%)")
    print("--------------------------------------------------")
    print(f"Avg Accel Duration Error            : {mean_acc_err:.2f}%")
    print(f"Avg Decel Duration Error            : {mean_dec_err:.2f}%")
    print(f"Avg Total Stroke Duration Error     : {mean_tot_err:.2f}%")
    print("==================================================")

if __name__ == "__main__":
    main()
