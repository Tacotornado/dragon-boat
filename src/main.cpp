#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <numeric>
#include "Filter.hpp"
#include "StrokeSegmenter.hpp"

// Simple CSV row parser to handle possible quotes and trailing whitespaces
std::vector<std::string> parse_csv_line(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ',')) {
        result.push_back(cell);
    }
    return result;
}

int main() {
    std::string raw_csv_path = "data/Canoe_xsens_dot_raw.csv";
    std::string out_csv_path = "reference/stroke_metrics_cpp.csv";

    std::ifstream file(raw_csv_path);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open raw data file: " << raw_csv_path << std::endl;
        return 1;
    }

    std::string header_line;
    if (!std::getline(file, header_line)) {
        std::cerr << "Error: Empty raw data file: " << raw_csv_path << std::endl;
        return 1;
    }

    // Identify columns index based on header
    std::vector<std::string> headers = parse_csv_line(header_line);
    int idx_time = -1, idx_ax = -1, idx_ay = -1, idx_az = -1;
    int idx_gx = -1, idx_gy = -1, idx_gz = -1;

    for (size_t i = 0; i < headers.size(); ++i) {
        std::string h = headers[i];
        if (h == "SampleTimeFine") idx_time = i;
        else if (h == "Acc_X") idx_ax = i;
        else if (h == "Acc_Y") idx_ay = i;
        else if (h == "Acc_Z") idx_az = i;
        else if (h == "Gyr_X") idx_gx = i;
        else if (h == "Gyr_Y") idx_gy = i;
        else if (h == "Gyr_Z") idx_gz = i;
    }

    if (idx_time == -1 || idx_ax == -1 || idx_ay == -1 || idx_az == -1 ||
        idx_gx == -1 || idx_gy == -1 || idx_gz == -1) {
        std::cerr << "Error: Missing required columns in raw CSV header." << std::endl;
        return 1;
    }

    // Instantiate processing pipeline
    const double FS = 120.0;
    const double DEFAULT_CUTOFF = 10.0;
    const int FILTER_ORDER = 4;

    Canoe::CausalAdaptiveFilter filter_ax(FILTER_ORDER, FS, DEFAULT_CUTOFF);
    Canoe::CausalAdaptiveFilter filter_ay(FILTER_ORDER, FS, DEFAULT_CUTOFF);
    Canoe::CausalAdaptiveFilter filter_az(FILTER_ORDER, FS, DEFAULT_CUTOFF);
    Canoe::CausalAdaptiveFilter filter_gx(FILTER_ORDER, FS, DEFAULT_CUTOFF);
    Canoe::CausalAdaptiveFilter filter_gy(FILTER_ORDER, FS, DEFAULT_CUTOFF);
    Canoe::CausalAdaptiveFilter filter_gz(FILTER_ORDER, FS, DEFAULT_CUTOFF);

    Canoe::PsdEstimator psd_estimator(360, 512, FS);
    Canoe::StrokeSegmenter segmenter(0.15, 0.08, FS);

    std::vector<Canoe::StrokeResult> detected_strokes;
    std::string line;
    int total_samples = 0;
    double current_cutoff = DEFAULT_CUTOFF;

    std::string filt_csv_path = "reference/filtered_imu_data.csv";
    std::ofstream filt_file(filt_csv_path);
    if (!filt_file.is_open()) {
        std::cerr << "Error: Could not open output file: " << filt_csv_path << std::endl;
        return 1;
    }
    filt_file << "SampleTimeFine,Acc_X,Acc_Y,Acc_Z,Gyr_X,Gyr_Y,Gyr_Z,Cutoff\n";

    std::cout << "Running C++ IMU pipeline..." << std::endl;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = parse_csv_line(line);
        if (row.size() < 8) continue;

        try {
            uint64_t timestamp_us = std::stoull(row[idx_time]);
            double ax_raw = std::stod(row[idx_ax]);
            double ay_raw = std::stod(row[idx_ay]);
            double az_raw = std::stod(row[idx_az]);
            double gx_raw = std::stod(row[idx_gx]);
            double gy_raw = std::stod(row[idx_gy]);
            double gz_raw = std::stod(row[idx_gz]);

            // 1. Update PSD & adapt cutoff frequency
            double cutoff = DEFAULT_CUTOFF;
            bool cutoff_updated = psd_estimator.process_sample(ay_raw, total_samples, cutoff);
            if (cutoff_updated) {
                current_cutoff = cutoff;
                filter_ax.update_cutoff(cutoff);
                filter_ay.update_cutoff(cutoff);
                filter_az.update_cutoff(cutoff);
                filter_gx.update_cutoff(cutoff);
                filter_gy.update_cutoff(cutoff);
                filter_gz.update_cutoff(cutoff);
            }

            // 2. Filter sample-by-sample causally
            Canoe::ImuSample filtered;
            filtered.timestamp_us = timestamp_us;
            filtered.sample_idx = total_samples;
            filtered.acc_x = filter_ax.process_sample(ax_raw);
            filtered.acc_y = filter_ay.process_sample(ay_raw);
            filtered.acc_z = filter_az.process_sample(az_raw);
            filtered.gyr_x = filter_gx.process_sample(gx_raw);
            filtered.gyr_y = filter_gy.process_sample(gy_raw);
            filtered.gyr_z = filter_gz.process_sample(gz_raw);

            // Write sample-by-sample filtered data
            filt_file << std::fixed << std::setprecision(16);
            filt_file << timestamp_us << ","
                      << filtered.acc_x << ","
                      << filtered.acc_y << ","
                      << filtered.acc_z << ","
                      << filtered.gyr_x << ","
                      << filtered.gyr_y << ","
                      << filtered.gyr_z << ","
                      << current_cutoff << "\n";

            // 3. Segment and calculate metrics
            auto opt_stroke = segmenter.process_sample(filtered);
            if (opt_stroke.has_value()) {
                detected_strokes.push_back(opt_stroke.value());
            }

            total_samples++;
        } catch (const std::exception& e) {
            // Skip invalid rows silently
        }
    }

    file.close();
    filt_file.close();

    std::cout << "Finished processing. Total samples: " << total_samples << std::endl;
    std::cout << "Detected " << detected_strokes.size() << " strokes." << std::endl;

    // Save metrics to output CSV
    std::ofstream out_file(out_csv_path);
    if (!out_file.is_open()) {
        std::cerr << "Error: Could not open output file: " << out_csv_path << std::endl;
        return 1;
    }

    // Write header
    out_file << "stroke_id,accel_duration_s,accY_pos_peak,t_to_pos_peak_pct,RAD,accel_velocity_integral,"
             << "accel_MAC_accX,accel_MAC_gyrX,accel_MAC_gyrY,accel_MAC_gyrZ,"
             << "decel_duration_s,accY_neg_peak,t_to_neg_peak_pct,decel_velocity_integral,"
             << "decel_MAC_accX,decel_MAC_gyrX,decel_MAC_gyrY,decel_MAC_gyrZ,"
             << "stroke_gyrX_pos_peak,stroke_gyrX_neg_peak,stroke_MAC_accX,"
             << "stroke_MAC_gyrX,stroke_MAC_gyrY,stroke_MAC_gyrZ,"
             << "TotalTime_s,TimeToPosPeak_s,TimeToNegPeak_s,Area_Positive,Area_Negative,"
             << "accel_start_idx,accel_end_idx,decel_start_idx,decel_end_idx,pos_peak_idx,neg_peak_idx\n";

    // Write rows with full precision to match Python
    out_file << std::fixed << std::setprecision(16);
    for (const auto& s : detected_strokes) {
        out_file << s.stroke_id << ","
                 << s.accel.duration_s << ","
                 << s.accel.peak << ","
                 << s.accel.t_to_peak_pct << ","
                 << s.rad << ","
                 << s.accel.velocity_integral << ","
                 << s.accel.mac_acc_x << ","
                 << s.accel.mac_gyr_x << ","
                 << s.accel.mac_gyr_y << ","
                 << s.accel.mac_gyr_z << ","
                 << s.decel.duration_s << ","
                 << s.decel.peak << ","
                 << s.decel.t_to_peak_pct << ","
                 << s.decel.velocity_integral << ","
                 << s.decel.mac_acc_x << ","
                 << s.decel.mac_gyr_x << ","
                 << s.decel.mac_gyr_y << ","
                 << s.decel.mac_gyr_z << ","
                 << s.gyr_x_pos_peak << ","
                 << s.gyr_x_neg_peak << ","
                 << s.stroke_mac_acc_x << ","
                 << s.stroke_mac_gyr_x << ","
                 << s.stroke_mac_gyr_y << ","
                 << s.stroke_mac_gyr_z << ","
                 << s.total_time_s << ","
                 << s.time_to_pos_peak_s << ","
                 << s.time_to_neg_peak_s << ","
                 << s.area_positive << ","
                 << s.area_negative << ","
                 << s.accel_start_idx << ","
                 << s.accel_end_idx << ","
                 << s.decel_start_idx << ","
                 << s.decel_end_idx << ","
                 << s.pos_peak_idx << ","
                 << s.neg_peak_idx << "\n";
    }

    out_file.close();
    std::cout << "Saved C++ metrics output to: " << out_csv_path << std::endl;

    // Calculate and display statistics for a few core parameters to check parity
    if (detected_strokes.empty()) return 0;

    double sum_tot_time = 0.0, sum_pos_peak = 0.0, sum_neg_peak = 0.0, sum_rad = 0.0;
    double min_tot_time = 1e9, min_pos_peak = 1e9, min_neg_peak = 1e9, min_rad = 1e9;
    double max_tot_time = -1e9, max_pos_peak = -1e9, max_neg_peak = -1e9, max_rad = -1e9;

    for (const auto& s : detected_strokes) {
        sum_tot_time += s.total_time_s;
        sum_pos_peak += s.accel.peak;
        sum_neg_peak += s.decel.peak;
        sum_rad += s.rad;

        min_tot_time = std::min(min_tot_time, s.total_time_s);
        min_pos_peak = std::min(min_pos_peak, s.accel.peak);
        min_neg_peak = std::min(min_neg_peak, s.decel.peak);
        min_rad = std::min(min_rad, s.rad);

        max_tot_time = std::max(max_tot_time, s.total_time_s);
        max_pos_peak = std::max(max_pos_peak, s.accel.peak);
        max_neg_peak = std::max(max_neg_peak, s.decel.peak);
        max_rad = std::max(max_rad, s.rad);
    }

    size_t sz = detected_strokes.size();
    std::cout << "\nC++ Pipeline Summary Statistics:" << std::endl;
    std::cout << std::setprecision(4);
    std::cout << "TotalTime (s): Mean=" << (sum_tot_time/sz) << ", Min=" << min_tot_time << ", Max=" << max_tot_time << std::endl;
    std::cout << "Pos Peak (m/s²): Mean=" << (sum_pos_peak/sz) << ", Min=" << min_pos_peak << ", Max=" << max_pos_peak << std::endl;
    std::cout << "Neg Peak (m/s²): Mean=" << (sum_neg_peak/sz) << ", Min=" << min_neg_peak << ", Max=" << max_neg_peak << std::endl;
    std::cout << "RAD: Mean=" << (sum_rad/sz) << ", Min=" << min_rad << ", Max=" << max_rad << std::endl;

    return 0;
}
