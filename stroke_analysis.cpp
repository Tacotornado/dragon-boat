//stroke_analysis.cpp
#include <bits/stdc++.h>
#include <filesystem>
#include "matplotlibcpp.h"
#include "filters.h"
#include "signal_features.h"
#include "stroke_metrics.h"
#include "math_utils.h"
#include <xlnt/xlnt.hpp>

namespace plt = matplotlibcpp;
namespace fs = std::filesystem;

// ---- StrokeCycle struct ----
struct StrokeCycle {
    int StrokeID;
    int StartZero, PosPeak, MidZero, NegPeak, EndZero;
    double TotalTime, TimeToPosPeak, TimeToNegPeak;
    double Area_Positive, Area_Negative;
};

struct AccPhaseMetrics {
    int StrokeID;
    double Duration;
    double AccY_PosPeak;
    double TimeToPosPeakPercent;
    double RAD;
    double VelocityChange;
    double MeanAbsAccX;
    double MeanAbsGyrX;
    double MeanAbsGyrY;
    double MeanAbsGyrZ;
};

struct DecPhaseMetrics {
    int StrokeID;
    double Duration;
    double AccY_NegPeak;
    double TimeToNegPeakPercent;
    double VelocityChange;
    double MeanAbsAccX;
    double MeanAbsGyrX;
    double MeanAbsGyrY;
    double MeanAbsGyrZ;
};

struct StrokePhaseMetrics {
    int StrokeID;
    double GyrX_PosPeak;
    double GyrX_NegPeak;
    double MeanAbsAccX;
    double MeanAbsGyrX;
    double MeanAbsGyrY;
    double MeanAbsGyrZ;
};

int main() {
    try {
        std::string RESULT_DIR = "results_full_strokes_A";
        std::string STROKE_PLOT_DIR = RESULT_DIR + "/stroke_plots";
        fs::create_directories(STROKE_PLOT_DIR);

        // --- Load Data from CSV ---
        std::ifstream fin("data/Canoe_xsens_dot.csv");
        if (!fin.is_open()) {
            std::cerr << "Failed to open CSV file!\n";
            return 1;
        }

        std::string line;
        std::getline(fin, line);
        std::stringstream ss(line);
        std::string col;
        std::vector<std::string> headers;
        while (std::getline(ss, col, ',')) headers.push_back(col);

        // --- Find column indices ---
        auto find_idx = [&](const std::string &name) {
            for (size_t i = 0; i < headers.size(); ++i)
                if (headers[i] == name) return (int)i;
            return -1;
        };

        int idx_AccX = find_idx("Acc_X");
        int idx_AccY = find_idx("Acc_Y");
        int idx_GyrX = find_idx("Gyr_X");
        int idx_GyrY = find_idx("Gyr_Y");
        int idx_GyrZ = find_idx("Gyr_Z");

        if (idx_AccY < 0) {
            std::cerr << "Acc_Y column not found!\n";
            return 1;
        }

        // --- Declare all channel vectors ---
        std::vector<double> accX, acc_y, gyrX, gyrY, gyrZ;


        // Read CSV rows and load full columns
        while (std::getline(fin, line)) {
            std::stringstream row(line);
            std::string cell;
            int col_index = 0;
            std::vector<std::string> row_cells;
            while (std::getline(row, cell, ',')) {
                row_cells.push_back(cell);
            }
            // guard number of columns
            if ((int)row_cells.size() <= std::max({idx_AccX, idx_AccY, idx_GyrX, idx_GyrY, idx_GyrZ})) {
                accX.push_back(0.0);
                acc_y.push_back(0.0);
                gyrX.push_back(0.0);
                gyrY.push_back(0.0);
                gyrZ.push_back(0.0);
                continue;
            }
            auto toD = [&](int idx){
                try { return std::stod(row_cells[idx]); }
                catch(...) { return 0.0; }
            };
            accX.push_back(toD(idx_AccX));
            acc_y.push_back(toD(idx_AccY));
            gyrX.push_back(toD(idx_GyrX));
            gyrY.push_back(toD(idx_GyrY));
            gyrZ.push_back(toD(idx_GyrZ));
        }

            // --- Filtering ---
            double fs_rate = 120.0; // Hz
            double cutoff = 5.0;   // Hz
            ButterworthLowpass lp_filter(cutoff, fs_rate, 2);
            auto acc_y_smooth = lp_filter.apply(acc_y);

            // --- Peaks ---
            auto pos_peaks = find_peaks_prominence(acc_y_smooth, 5, 0.02);
            auto neg_peaks = find_peaks_prominence_neg(acc_y_smooth, 5, 0.02);
            

            std::cout << "Pos peaks: " << pos_peaks.size() << " Neg peaks: " << neg_peaks.size() << "\n";

            // --- Zero crossings (fractional positions) ---
            std::vector<double> all_zeros = find_all_zero_crossings(acc_y_smooth, 0.1, 7, 0.01, (int)(0.15*fs_rate), 0.01);
            std::cout << "Zero crossings: " << all_zeros.size() << "\n";

            // --- Stroke cycles ---
            std::vector<StrokeCycle> stroke_cycles;
            std::vector<AccPhaseMetrics> acc_metrics;
            std::vector<DecPhaseMetrics> dec_metrics;
            std::vector<StrokePhaseMetrics> stroke_phase_metrics;

// --- Stroke cycles (using fractional zeros for timing, integer indices for slicing) ---
for (size_t i = 0; i + 2 < all_zeros.size(); ++i) {
    double z1d = all_zeros[i];
    double z2d = all_zeros[i+1];
    double z3d = all_zeros[i+2];

    // integer bounds for slicing / peak search 
    int z1 = std::max(0, (int)std::floor(z1d));
    int z2 = std::min((int)acc_y_smooth.size()-1, (int)std::round(z2d));
    int z3 = std::min((int)acc_y_smooth.size()-1, (int)std::ceil(z3d));

    if (z1 >= z2 || z2 >= z3) continue;

    // Find positive peaks inside (z1, z2)
    std::vector<double> pos_in_first;
    for (double p : pos_peaks)
        if (p > z1 && p < z2)
            pos_in_first.push_back((int)std::round(p));

    // Find negative peaks inside (z2, z3)
    std::vector<double> neg_in_second;
    for (double p : neg_peaks)
        if (p > z2 && p < z3)
            neg_in_second.push_back((int)std::round(p));


    if (pos_in_first.empty() || neg_in_second.empty()) continue;

    // strongest positive in first half
    int pos_idx = *std::max_element(pos_in_first.begin(), pos_in_first.end(),
        [&](int a,int b){ return acc_y_smooth[a] < acc_y_smooth[b]; });

    // strongest negative in second half
    int neg_idx = *std::min_element(neg_in_second.begin(), neg_in_second.end(),
        [&](int a,int b){ return acc_y_smooth[a] < acc_y_smooth[b]; });

    if (acc_y_smooth[pos_idx] <= 0 || acc_y_smooth[neg_idx] >= 0) continue;

    // --- Timing using fractional zero positions (higher precision) ---
    double totalTime = (z3d - z1d) / fs_rate;
    double timeToPosPeak = (pos_idx - z1d) / fs_rate;
    double timeToNegPeak = (neg_idx - z2d) / fs_rate;

    // Compute areas using fractional zero crossings 
    auto slice_with_fractional = [&](double start, double end) {
        std::vector<double> x, y;
        int i_start = (int)std::floor(start);
        int i_end   = (int)std::ceil(end);

        for (int k = i_start; k <= i_end; ++k) {
            double xi = k;
            double yi = interp_at(acc_y_smooth, xi);  // fractional y interpolation
            x.push_back(xi);
            y.push_back(yi);
        }
        return std::make_pair(x, y);
    };

    // Positive area (z1d → z2d)
    auto [x_pos, y_pos] = slice_with_fractional(z1d, z2d);
    for (double &v : y_pos) v = std::max(0.0, v);
    double pos_area = trapz(y_pos, 1.0/fs_rate);

    // Negative area (z2d → z3d)
    auto [x_neg, y_neg] = slice_with_fractional(z2d, z3d);
    for (double &v : y_neg) v = std::min(0.0, v);
    double neg_area_mag = std::abs(trapz(y_neg, 1.0/fs_rate));

    // assign StrokeID and use integer indices for Start/Mid/End to preserve current downstream logic
    StrokeCycle sc{
        (int)stroke_cycles.size() + 1,
        z1, pos_idx, z2, neg_idx, z3,
        totalTime, timeToPosPeak, timeToNegPeak,
        pos_area, neg_area_mag
    };

    stroke_cycles.push_back(sc);

    // --- compute metrics as before; use integer z1,z2,z3 for slicing
    const StrokeCycle &last = stroke_cycles.back();
    int sid = last.StrokeID;

    // Acceleration phase
    double dur_acc = (z2 - z1) / fs_rate;
    double accY_peak = acc_y_smooth[last.PosPeak];
    double time_to_peak = (last.PosPeak - z1) / fs_rate;
    double percent_time = (time_to_peak / std::max(1e-12, dur_acc)) * 100.0;
    double RAD = accY_peak / std::max(1e-12, time_to_peak);
    double vel_change_acc = trapz(std::vector<double>(acc_y_smooth.begin()+z1, acc_y_smooth.begin()+z2), 1.0/fs_rate);

    double mean_abs_accx_acc = mean_abs_diff(std::vector<double>(accX.begin()+z1, accX.begin()+z2));
    double mean_abs_gyx_acc  = mean_abs_diff(std::vector<double>(gyrX.begin()+z1, gyrX.begin()+z2));
    double mean_abs_gyy_acc  = mean_abs_diff(std::vector<double>(gyrY.begin()+z1, gyrY.begin()+z2));
    double mean_abs_gyz_acc  = mean_abs_diff(std::vector<double>(gyrZ.begin()+z1, gyrZ.begin()+z2));

    acc_metrics.push_back({
        sid, dur_acc, accY_peak, percent_time, RAD, vel_change_acc,
        mean_abs_accx_acc, mean_abs_gyx_acc, mean_abs_gyy_acc, mean_abs_gyz_acc
    });

    // Deceleration phase
    double dur_dec = (z3 - z2) / fs_rate;
    double accY_negpeak = acc_y_smooth[last.NegPeak];
    double time_to_neg = (last.NegPeak - z2) / fs_rate;
    double percent_neg_time = (time_to_neg / std::max(1e-12, dur_dec)) * 100.0;
    double vel_change_dec = trapz(std::vector<double>(acc_y_smooth.begin()+z2, acc_y_smooth.begin()+z3), 1.0/fs_rate);

    double mean_abs_accx_dec = mean_abs_diff(std::vector<double>(accX.begin()+z2, accX.begin()+z3));
    double mean_abs_gyx_dec  = mean_abs_diff(std::vector<double>(gyrX.begin()+z2, gyrX.begin()+z3));
    double mean_abs_gyy_dec  = mean_abs_diff(std::vector<double>(gyrY.begin()+z2, gyrY.begin()+z3));
    double mean_abs_gyz_dec  = mean_abs_diff(std::vector<double>(gyrZ.begin()+z2, gyrZ.begin()+z3));

    dec_metrics.push_back({
        sid, dur_dec, accY_negpeak, percent_neg_time, vel_change_dec,
        mean_abs_accx_dec, mean_abs_gyx_dec, mean_abs_gyy_dec, mean_abs_gyz_dec
    });

    // Stroke phase metrics 
    auto gyrx_slice = std::vector<double>(gyrX.begin()+z1, gyrX.begin()+z3);
    double gyrx_pos_peak = *std::max_element(gyrx_slice.begin(), gyrx_slice.end());
    double gyrx_neg_peak = *std::min_element(gyrx_slice.begin(), gyrx_slice.end());
    double mean_abs_accx_stroke = mean_abs_diff(std::vector<double>(accX.begin()+z1, accX.begin()+z3));
    double mean_abs_gyx_stroke  = mean_abs_diff(std::vector<double>(gyrX.begin()+z1, gyrX.begin()+z3));
    double mean_abs_gyy_stroke  = mean_abs_diff(std::vector<double>(gyrY.begin()+z1, gyrY.begin()+z3));
    double mean_abs_gyz_stroke  = mean_abs_diff(std::vector<double>(gyrZ.begin()+z1, gyrZ.begin()+z3));

    stroke_phase_metrics.push_back({
    sid, gyrx_pos_peak, gyrx_neg_peak,
    mean_abs_accx_stroke, mean_abs_gyx_stroke,
    mean_abs_gyy_stroke, mean_abs_gyz_stroke
});
} 


        std::cout << "Saved " << stroke_cycles.size() << " stroke cycles\n";
        
        xlnt::workbook wb;

        // --- Sheet 1: Acceleration Phase ---
        auto ws1 = wb.active_sheet();
        ws1.title("Acceleration Phase");

        std::vector<std::string> headers1 = {
            "StrokeID","Duration","ACC-Y+Peak","TimeToPeak(%)","RAD",
            "ΔVelocity","|ACC-X|","|Gyr-X|","|Gyr-Y|","|Gyr-Z|"
        };
        for (size_t i = 0; i < headers1.size(); ++i)
            ws1.cell(1, i + 1).value(headers1[i]);

        int row = 2;
        for (auto &m : acc_metrics) {
            ws1.cell(row, 1).value(m.StrokeID);
            ws1.cell(row, 2).value(m.Duration);
            ws1.cell(row, 3).value(m.AccY_PosPeak);
            ws1.cell(row, 4).value(m.TimeToPosPeakPercent);
            ws1.cell(row, 5).value(m.RAD);
            ws1.cell(row, 6).value(m.VelocityChange);
            ws1.cell(row, 7).value(m.MeanAbsAccX);
            ws1.cell(row, 8).value(m.MeanAbsGyrX);
            ws1.cell(row, 9).value(m.MeanAbsGyrY);
            ws1.cell(row, 10).value(m.MeanAbsGyrZ);
            ++row;
        }

        // --- Sheet 2: Deceleration Phase ---
        auto ws2 = wb.create_sheet();
        ws2.title("Deceleration Phase");
        std::vector<std::string> headers2 = {
            "StrokeID","Duration","ACC-Y−Peak","TimeToNegPeak(%)",
            "ΔVelocity","|ACC-X|","|Gyr-X|","|Gyr-Y|","|Gyr-Z|"
        };
        for (size_t i = 0; i < headers2.size(); ++i)
            ws2.cell(1, i + 1).value(headers2[i]);

        int r2 = 2;
        for (auto &m : dec_metrics) {
            ws2.cell(r2, 1).value(m.StrokeID);
            ws2.cell(r2, 2).value(m.Duration);
            ws2.cell(r2, 3).value(m.AccY_NegPeak);
            ws2.cell(r2, 4).value(m.TimeToNegPeakPercent);
            ws2.cell(r2, 5).value(m.VelocityChange);
            ws2.cell(r2, 6).value(m.MeanAbsAccX);
            ws2.cell(r2, 7).value(m.MeanAbsGyrX);
            ws2.cell(r2, 8).value(m.MeanAbsGyrY);
            ws2.cell(r2, 9).value(m.MeanAbsGyrZ);
            ++r2;
        }


        // --- Sheet 3: Stroke Phase ---
        auto ws3 = wb.create_sheet();
        ws3.title("Stroke Phase");
        std::vector<std::string> headers3 = {
            "StrokeID","GyrX+Peak","GyrX−Peak","|ACC-X|","|Gyr-X|","|Gyr-Y|","|Gyr-Z|"
        };
        for (size_t i = 0; i < headers3.size(); ++i)
            ws3.cell(1, i + 1).value(headers3[i]);

        int r3 = 2;
        for (auto &m : stroke_phase_metrics) {
            ws3.cell(r3, 1).value(m.StrokeID);
            ws3.cell(r3, 2).value(m.GyrX_PosPeak);
            ws3.cell(r3, 3).value(m.GyrX_NegPeak);
            ws3.cell(r3, 4).value(m.MeanAbsAccX);
            ws3.cell(r3, 5).value(m.MeanAbsGyrX);
            ws3.cell(r3, 6).value(m.MeanAbsGyrY);
            ws3.cell(r3, 7).value(m.MeanAbsGyrZ);
            ++r3;
        }


        wb.save("results_full_strokes_A/stroke_phase_metrics.xlsx");

        // --- Individual stroke plots (start & end exactly at zero crossings) ---
        for (size_t i = 0; i < stroke_cycles.size(); ++i) {
            auto &sc = stroke_cycles[i];

            // Get fractional zero crossings for this stroke
            double z1d = all_zeros[i * 2];       // start zero (fractional)
            double z2d = all_zeros[i * 2 + 1];   // mid zero
            double z3d = all_zeros[i * 2 + 2];   // end zero (fractional)

            // Create fine-grained sample grid (e.g., 0.1 sample steps)
            double step = 0.1;
            std::vector<double> x, y;
            for (double xi = z1d; xi <= z3d; xi += step) {
                x.push_back(xi);
                y.push_back(interp_at(acc_y_smooth, xi)); // linear interpolation
            }

            plt::figure_size(800, 400);
            plt::plot(x, y, "b-");

            // Mark fractional zero crossings
            plt::plot({z1d}, {interp_at(acc_y_smooth, z1d)}, "kx"); // start
            plt::plot({z2d}, {interp_at(acc_y_smooth, z2d)}, "kx"); // mid
            plt::plot({z3d}, {interp_at(acc_y_smooth, z3d)}, "kx"); // end

            // Mark peaks (still integer indices)
            plt::plot({(double)sc.PosPeak}, {acc_y_smooth[sc.PosPeak]}, "ro");
            plt::plot({(double)sc.NegPeak}, {acc_y_smooth[sc.NegPeak]}, "bv");

            plt::title("Stroke " + std::to_string(sc.StrokeID) +
                    " [zero: " + std::to_string(z1d) + " → " + std::to_string(z3d) + "]");
            plt::xlabel("Sample index");
            plt::ylabel("Acc_Y");

            std::string fname = STROKE_PLOT_DIR + "/stroke_" + std::to_string(sc.StrokeID) + ".png";
            plt::save(fname);
            plt::close();
        }

        // --- Save CSV ---
        fs::create_directories(RESULT_DIR);
        std::ofstream out(RESULT_DIR+"/stroke_full_cycles.csv");
        out << "StrokeID,StartZero,PosPeak,MidZero,NegPeak,EndZero,TotalTime,TimeToPosPeak,TimeToNegPeak,Area_Positive,Area_Negative\n";
        for (auto &sc: stroke_cycles) {
            out << sc.StrokeID << "," << sc.StartZero << "," << sc.PosPeak << "," 
                << sc.MidZero << "," << sc.NegPeak << "," << sc.EndZero << ","
                << sc.TotalTime << "," << sc.TimeToPosPeak << "," << sc.TimeToNegPeak << ","
                << sc.Area_Positive << "," << sc.Area_Negative << "\n";
        }
        out.close();

        // --- Plot overview ---
        plt::figure_size(1600,600);
        plt::plot(acc_y_smooth);

        std::vector<double> pos_x, pos_y, neg_x, neg_y;
        for (auto &sc : stroke_cycles) {
            pos_x.push_back((double)sc.PosPeak);
            pos_y.push_back(acc_y_smooth[sc.PosPeak]);
            neg_x.push_back((double)sc.NegPeak);
            neg_y.push_back(acc_y_smooth[sc.NegPeak]);
        }
        if (!pos_x.empty()) plt::plot(pos_x, pos_y, "ro");
        if (!neg_x.empty()) plt::plot(neg_x, neg_y, "bv");

        // --- Zero crossing markers (plot at exact interpolated y; here we use interp_at to get y) ---
        std::vector<double> zero_x, zero_y;
        for (double zpos : all_zeros) {
            zero_x.push_back(zpos);
            // use interp_at to find fractional-sample y value
            zero_y.push_back(interp_at(acc_y_smooth, zpos));
        }
        if (!zero_x.empty()) plt::plot(zero_x, zero_y, "kx");

        plt::title("Detected stroke cycles");
        plt::xlabel("Sample index");
        plt::ylabel("Acc_Y");

        plt::save(RESULT_DIR+"/strokes_overview.png");
        plt::close();

    } catch(std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
}