#include <bits/stdc++.h>
#include <filesystem>
#include "matplotlibcpp.h"
#include "filters.h"

namespace plt = matplotlibcpp;
namespace fs = std::filesystem;


// ---- Utility functions ----
double mean(const std::vector<double> &v) {
    if (v.empty()) return 0.0;
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

double mean_abs_diff(const std::vector<double> &v) {
    if (v.size() < 2) return 0.0;
    double sum = 0.0;
    for (size_t i = 1; i < v.size(); i++) sum += std::abs(v[i] - v[i-1]);
    return sum / (v.size()-1);
}

double trapz(const std::vector<double> &y, double dx) {
    if (y.size() < 2) return 0.0;
    double s = 0.0;
    for (size_t i = 1; i < y.size(); i++)
        s += 0.5 * (y[i-1] + y[i]) * dx;
    return s;
}

// ---- Peak detection with prominence ----
std::vector<int> find_peaks_prominence(const std::vector<double> &sig, int distance, double prominence) {
    std::vector<int> peaks;
    for (int i = distance; i < (int)sig.size()-distance; i++) {
        if (sig[i] > sig[i-1] && sig[i] > sig[i+1]) {
            bool ok = true;
            for (int j = 1; j <= distance; j++) {
                if (sig[i] < sig[i-j] || sig[i] < sig[i+j]) { ok = false; break; }
            }
            if (!ok) continue;

            double left_min = sig[i];
            for (int j = i; j >= 0; j--) left_min = std::min(left_min, sig[j]);
            double right_min = sig[i];
            for (int j = i; j < (int)sig.size(); j++) right_min = std::min(right_min, sig[j]);
            double prom = sig[i] - std::max(left_min, right_min);
            if (prom >= prominence) peaks.push_back(i);
        }
    }
    return peaks;
}

std::vector<int> find_peaks_prominence_neg(const std::vector<double> &sig, int distance, double prominence) {
    std::vector<int> peaks;
    for (int i = distance; i < (int)sig.size()-distance; i++) {
        if (sig[i] < sig[i-1] && sig[i] < sig[i+1]) {
            bool ok = true;
            for (int j = 1; j <= distance; j++) {
                if (sig[i] > sig[i-j] || sig[i] > sig[i+j]) { ok = false; break; }
            }
            if (!ok) continue;

            double left_max = sig[i];
            for (int j = i; j >= 0; j--) left_max = std::max(left_max, sig[j]);
            double right_max = sig[i];
            for (int j = i; j < (int)sig.size(); j++) right_max = std::max(right_max, sig[j]);
            double prom = std::min(left_max, right_max) - sig[i];
            if (prom >= prominence) peaks.push_back(i);
        }
    }
    return peaks;
}

// ---- Zero crossing detection ----
std::vector<int> find_all_zero_crossings(const std::vector<double> &sig,
                                         double min_amp, int window,
                                         double slope_thr, int min_spacing,
                                         double hysteresis) {
    std::vector<int> zeros;
    int last_z = -min_spacing;
    auto sign = [&](double x){ return (x > hysteresis) - (x < -hysteresis); };

    for (int i = 1; i < (int)sig.size(); i++) {
        int s_prev = sign(sig[i-1]);
        int s_curr = sign(sig[i]);
        if (s_prev == 0) s_prev = sign(sig[i-2 >= 0 ? i-2 : 0]); // stabilize near flat zero
        if (s_prev != 0 && s_curr != 0 && s_prev != s_curr) {
            double slope = sig[i] - sig[i-1];
            int lo = std::max(0, i-window);
            int hi = std::min((int)sig.size(), i+window);
            double local_max = *std::max_element(sig.begin()+lo, sig.begin()+hi);
            double local_min = *std::min_element(sig.begin()+lo, sig.begin()+hi);
            if ((local_max - local_min) > min_amp && std::abs(slope) > slope_thr) {
                if (i - last_z >= min_spacing) {
                    zeros.push_back(i);
                    last_z = i;
                }
            }
        }
    }
    return zeros;
}


// ---- StrokeCycle struct ----
struct StrokeCycle {
    int StrokeID;
    int StartZero, PosPeak, MidZero, NegPeak, EndZero;
    double TotalTime, TimeToPosPeak, TimeToNegPeak;
    double Area_Positive, Area_Negative;
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
        std::getline(fin, line); // header
        std::stringstream ss(line);
        std::string col;
        std::vector<std::string> headers;
        while (std::getline(ss, col, ',')) headers.push_back(col);

        // Find Acc_Y column
        int idx_AccY = -1;
        for (size_t i = 0; i < headers.size(); i++) {
            if (headers[i] == "Acc_Y") { idx_AccY = (int)i; break; }
        }
        if (idx_AccY < 0) {
            std::cerr << "Acc_Y column not found!\n";
            return 1;
        }

        // Read Acc_Y values
        std::vector<double> acc_y;
        while (std::getline(fin, line)) {
            std::stringstream row(line);
            std::string cell;
            int col_index = 0;
            while (std::getline(row, cell, ',')) {
                if (col_index == idx_AccY) {
                    try { acc_y.push_back(std::stod(cell)); }
                    catch (...) { acc_y.push_back(0.0); }
                    break;
                }
                col_index++;
            }
        }
        fin.close();
        std::cout << "Loaded samples: " << acc_y.size() << "\n";

        // --- Filtering ---
        double fs_rate = 120.0; // Hz
        double cutoff = 10.0;   // Hz
        ButterworthLowpass lp_filter(cutoff, fs_rate, 2);
        auto acc_y_smooth = lp_filter.apply(acc_y);

        // --- Peaks ---
        auto pos_peaks = find_peaks_prominence(acc_y_smooth, 5, 0.02);
        auto neg_peaks = find_peaks_prominence_neg(acc_y_smooth, 5, 0.02);
        

        std::cout << "Pos peaks: " << pos_peaks.size() << " Neg peaks: " << neg_peaks.size() << "\n";

        // --- Zero crossings ---
        auto all_zeros = find_all_zero_crossings(acc_y_smooth, 0.1, 7, 0.01, (int)(0.15*fs_rate), 0.01);
        std::cout << "Zero crossings: " << all_zeros.size() << "\n";

        // --- Stroke cycles ---
        std::vector<StrokeCycle> stroke_cycles;
        for (size_t i = 0; i + 2 < all_zeros.size(); i++) {   // fix guard: need z1,z2,z3
            int z1 = all_zeros[i];
            int z2 = all_zeros[i+1];
            int z3 = all_zeros[i+2];

            // Find peaks separately in each half
            std::vector<int> pos_in_first, neg_in_second;
            for (int p : pos_peaks) if (p > z1 && p < z2) pos_in_first.push_back(p);
            for (int p : neg_peaks) if (p > z2 && p < z3) neg_in_second.push_back(p);

            if (pos_in_first.empty() || neg_in_second.empty()) continue;

            // Strongest positive in first half
            int pos_idx = *std::max_element(pos_in_first.begin(), pos_in_first.end(),
                [&](int a,int b){ return acc_y_smooth[a] < acc_y_smooth[b]; });
            // Strongest negative in second half
            int neg_idx = *std::min_element(neg_in_second.begin(), neg_in_second.end(),
                [&](int a,int b){ return acc_y_smooth[a] < acc_y_smooth[b]; });

            // Sanity: signs
            if (acc_y_smooth[pos_idx] <= 0 || acc_y_smooth[neg_idx] >= 0) continue;

            double t_total = (z3 - z1)/fs_rate;
            double t_pos   = (pos_idx - z1)/fs_rate;
            double t_neg   = (neg_idx - z2)/fs_rate;

            // Areas: strict halves
            std::vector<double> segpos, segneg;
            segpos.reserve(z2 - z1 + 1);
            segneg.reserve(z3 - z2 + 1);
            for (int k = z1; k <= z2; k++) segpos.push_back(std::max(0.0, acc_y_smooth[k]));
            for (int k = z2; k <= z3; k++) segneg.push_back(std::min(0.0, acc_y_smooth[k]));
            double pos_area     = trapz(segpos, 1.0/fs_rate);
            double neg_area_mag = std::abs(trapz(segneg, 1.0/fs_rate));

            stroke_cycles.push_back({(int)stroke_cycles.size()+1, z1,pos_idx,z2,neg_idx,z3,
                                    t_total,t_pos,t_neg,pos_area,neg_area_mag});
        }

        std::cout << "Saved " << stroke_cycles.size() << " stroke cycles\n";
        
        // --- Individual stroke plots ---
        for (auto &sc : stroke_cycles) {
            // Extract segment
            std::vector<double> seg;
            for (int k = sc.StartZero; k <= sc.EndZero; k++)
                seg.push_back(acc_y_smooth[k]);

            // X-axis relative to segment start
            std::vector<double> x(seg.size());
            std::iota(x.begin(), x.end(), 0);

            plt::figure_size(800,400);
            plt::plot(x, seg, "b-");

            // Mark positive peak
            int pos_rel = sc.PosPeak - sc.StartZero;
            plt::plot({(double)pos_rel}, {acc_y_smooth[sc.PosPeak]}, "ro");

            // Mark mid zero
            int mid_rel = sc.MidZero - sc.StartZero;
            plt::plot({(double)mid_rel}, {0.0}, "kx");

            // Mark negative peak
            int neg_rel = sc.NegPeak - sc.StartZero;
            plt::plot({(double)neg_rel}, {acc_y_smooth[sc.NegPeak]}, "bv");

            // Mark start and end zeros
            plt::plot({0.0}, {0.0}, "kx");
            plt::plot({(double)(sc.EndZero - sc.StartZero)}, {0.0}, "kx");

            plt::title("Stroke " + std::to_string(sc.StrokeID));
            plt::xlabel("Sample offset");
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

        // --- Zero crossing markers ---
        std::vector<double> zero_x, zero_y;
        for (int z : all_zeros) {
            zero_x.push_back((double)z);
            zero_y.push_back(0.0);  // mark on x-axis
        }
        if (!zero_x.empty()) plt::plot(zero_x, zero_y, "kx"); // black 'x' marks

        plt::title("Detected stroke cycles");
        plt::xlabel("Sample index");
        plt::ylabel("Acc_Y");
        plt::save(RESULT_DIR+"/strokes_overview.png");
        plt::close();

    } catch(std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
}