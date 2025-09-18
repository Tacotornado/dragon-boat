// stroke_analysis.cpp
// Full line-by-line port of the provided Python script to C++.
// Dependencies:
//  - Eigen (matrix/vector math)
//  - xlnt (read/write .xlsx)
//  - matplotlib-cpp (plotting; requires Python + matplotlib installed)
//  - C++17 or later
// Build (example using g++):
//  g++ -std=c++17 stroke_analysis.cpp -I/path/to/eigen -I/path/to/xlnt/include -I/path/to/matplotlibcpp -lpython3.10 -o stroke_analysis
// Adjust python library version and include paths as needed.

#include "xlnt/xlnt.hpp"
#include <Eigen/Dense>
#include "matplotlibcpp.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <filesystem>
#include <sstream>

namespace plt = matplotlibcpp;
namespace fs = std::filesystem;
using Eigen::VectorXd;

// ---------- Utility helpers ----------
static double nan_val = std::numeric_limits<double>::quiet_NaN();

inline int clamp_idx(int i, int n) { return std::max(0, std::min(n-1, i)); }

std::vector<double> slice(const std::vector<double>& v, int a, int b) {
    // return v[a:b] inclusive of a and exclusive of b (like Python)
    std::vector<double> out;
    if (b <= a) return out;
    a = std::max(0, a);
    b = std::min((int)v.size(), b);
    out.insert(out.end(), v.begin()+a, v.begin()+b);
    return out;
}

// mean of vector
double mean_vec(const std::vector<double>& v) {
    if (v.empty()) return nan_val;
    double s = std::accumulate(v.begin(), v.end(), 0.0);
    return s / v.size();
}

// mean of absolute values of diffs
double mean_abs_diff(const std::vector<double>& v) {
    if (v.size() < 2) return 0.0;
    double s = 0.0;
    for (size_t i = 1; i < v.size(); ++i) s += std::abs(v[i] - v[i-1]);
    return s / (v.size()-1);
}

// trapz integration
double trapz(const std::vector<double>& y, double dx) {
    if (y.size() < 2) return 0.0;
    double area = 0.0;
    for (size_t i = 1; i < y.size(); ++i) {
        area += 0.5 * (y[i] + y[i-1]) * dx;
    }
    return area;
}

// ---------- Digital Butterworth design (bilinear transform) ----------
// Returns IIR coefficients b (numerator) and a (denominator) for lowpass butterworth
// using bilinear transform with pre-warping. This is a straightforward implementation
// for order up to small values. We compute using analog prototype poles then apply
// bilinear transform. For simplicity we use the cookbook-style cascaded biquad approach
// but here we produce 'a' and 'b' arrays for direct-form filtering using scipy-like design.

struct IIRCoeffs { std::vector<double> b, a; };

IIRCoeffs butter_lowpass_coeffs(double cutoff, double fs, int order=4) {
    // Normalize
    double nyq = 0.5 * fs;
    double normal_cutoff = cutoff / nyq; // 0..1
    if (normal_cutoff <= 0 || normal_cutoff >= 1) {
        throw std::runtime_error("cutoff must be between 0 and Nyquist");
    }

    // We'll use bilinear transform conversion with pre-warped analog frequency
    double warped = std::tan(M_PI * normal_cutoff) ; // prewarp
    // Create poles of normalized analog Butterworth
    std::vector<std::complex<double>> poles;
    for (int k = 0; k < order; ++k) {
        double theta = M_PI * (2.0*k + 1.0) / (2.0*order);
        std::complex<double> pk = std::polar(1.0, theta);
        pk = std::complex<double>(-std::real(pk), -std::imag(pk)); // left-half plane
        poles.push_back(pk);
    }

    // Convert analog prototype to lowpass with cutoff 'warped'
    for (auto &p : poles) p *= warped;

    // Bilinear transform: s = (2/T)*(z-1)/(z+1), T = 1/fs. Build digital filter from poles
    // We form transfer function H(z) = K * product (z - z_k) / product (z - p_zk)
    // For simplicity and small order, expand polynomial from poles and zeros.

    // zeros are at z = -1 (mapped from infinity) with multiplicity = order
    // compute polynomial from poles in z-domain: map each analog pole s_p to z-pole via
    // z = (1 + s_p * T/2) / (1 - s_p * T/2)
    double T = 1.0 / fs;
    std::vector<std::complex<double>> z_poles;
    for (auto &p : poles) {
        std::complex<double> z = (1.0 + p * T / 2.0) / (1.0 - p * T / 2.0);
        z_poles.push_back(z);
    }

    // Build denominator polynomial from z_poles
    std::vector<std::complex<double>> a_poly = {1.0};
    for (auto &zp : z_poles) {
        std::vector<std::complex<double>> next(a_poly.size()+1);
        for (size_t i = 0; i < a_poly.size(); ++i) {
            next[i] += a_poly[i] * (-zp);
            next[i+1] += a_poly[i];
        }
        a_poly.swap(next);
    }

    // Build numerator polynomial: zeros at z=-1 (multiplicity=order)
    std::vector<std::complex<double>> b_poly = {1.0};
    for (int i = 0; i < order; ++i) {
        std::vector<std::complex<double>> next(b_poly.size()+1);
        for (size_t j = 0; j < b_poly.size(); ++j) {
            next[j] += b_poly[j] * (1.0);   // multiply by (z - (-1)) = (z+1)
            next[j+1] += b_poly[j] * (1.0);
        }
        b_poly.swap(next);
    }

    // Gain normalization: set DC gain to 1
    // Evaluate numerator and denominator at z=1
    std::complex<double> num(0.0,0.0), den(0.0,0.0);
    for (size_t i = 0; i < b_poly.size(); ++i) num += b_poly[i] * std::pow(1.0, (int)(b_poly.size()-1-i));
    for (size_t i = 0; i < a_poly.size(); ++i) den += a_poly[i] * std::pow(1.0, (int)(a_poly.size()-1-i));
    std::complex<double> K = den / num;
    for (auto &c : b_poly) c *= K;

    // Convert to real coefficients (imag small)
    std::vector<double> b_real(b_poly.size()), a_real(a_poly.size());
    for (size_t i = 0; i < b_poly.size(); ++i) b_real[i] = b_poly[i].real();
    for (size_t i = 0; i < a_poly.size(); ++i) a_real[i] = a_poly[i].real();

    return {b_real, a_real};
}

// IIR filtering (direct form using a and b like scipy.signal.lfilter)
std::vector<double> lfilter(const std::vector<double>& b, const std::vector<double>& a, const std::vector<double>& x) {
    int n = x.size();
    int nb = b.size();
    int na = a.size();
    std::vector<double> y(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double acc = 0.0;
        for (int j = 0; j < nb; ++j) {
            if (i-j >= 0) acc += b[j] * x[i-j];
        }
        for (int j = 1; j < na; ++j) {
            if (i-j >= 0) acc -= a[j] * y[i-j];
        }
        acc /= a[0];
        y[i] = acc;
    }
    return y;
}

// filtfilt: forward-backward filtering to approximate zero-phase
std::vector<double> filtfilt(const std::vector<double>& b, const std::vector<double>& a, const std::vector<double>& x) {
    // Simple approach: pad by reflecting half of the signal, run lfilter forward then reverse and filter
    int n = x.size();
    int padlen = 3 * std::max((int)b.size(), (int)a.size());
    if (n <= padlen) padlen = n;
    std::vector<double> extended; extended.reserve(n + 2*padlen);
    // reflect start
    for (int i = padlen-1; i >= 0; --i) extended.push_back(x[i]);
    // original
    for (int i = 0; i < n; ++i) extended.push_back(x[i]);
    // reflect end
    for (int i = n-1; i >= n-padlen; --i) extended.push_back(x[i]);

    auto y = lfilter(b, a, extended);
    // reverse
    std::reverse(y.begin(), y.end());
    y = lfilter(b, a, y);
    std::reverse(y.begin(), y.end());

    // remove padding
    std::vector<double> out(y.begin()+padlen, y.begin()+padlen+n);
    return out;
}

std::vector<double> lowpass_filter(const std::vector<double>& data, double cutoff=10.0, double fs=120.0, int order=4) {
    auto coeffs = butter_lowpass_coeffs(cutoff, fs, order);
    return filtfilt(coeffs.b, coeffs.a, data);
}

// ---------- Peak detection with prominence & distance ----------

std::vector<int> find_peaks_prominence(const std::vector<double>& data, int distance=5, double prominence=0.8) {
    int n = data.size();
    std::vector<int> peaks;
    for (int i = 1; i < n-1; ++i) {
        if (data[i] > data[i-1] && data[i] > data[i+1]) peaks.push_back(i);
    }
    // enforce distance: keep highest peaks
    std::vector<char> keep(peaks.size(), 1);
    for (size_t i = 0; i < peaks.size(); ++i) {
        if (!keep[i]) continue;
        for (size_t j = i+1; j < peaks.size(); ++j) {
            if (!keep[j]) continue;
            if (std::abs(peaks[i] - peaks[j]) <= distance) {
                if (data[peaks[i]] >= data[peaks[j]]) keep[j] = 0; else keep[i] = 0;
            }
        }
    }
    std::vector<int> filtered;
    for (size_t i = 0; i < peaks.size(); ++i) if (keep[i]) filtered.push_back(peaks[i]);

    // compute prominence: for each peak, find highest minima on both sides until a higher peak is encountered
    std::vector<int> final_peaks;
    for (int p : filtered) {
        // left min
        double left_min = data[p];
        for (int i = p; i >= 0; --i) {
            left_min = std::min(left_min, data[i]);
            if (data[i] > data[p]) break;
        }
        double right_min = data[p];
        for (int i = p; i < n; ++i) {
            right_min = std::min(right_min, data[i]);
            if (data[i] > data[p]) break;
        }
        double baseline = std::max(left_min, right_min);
        double prom = data[p] - baseline;
        if (prom >= prominence) final_peaks.push_back(p);
    }
    return final_peaks;
}

std::vector<int> find_peaks_prominence_neg(const std::vector<double>& data, int distance=5, double prominence=0.8) {
    // find peaks of negative signal
    std::vector<double> neg(data.size());
    for (size_t i = 0; i < data.size(); ++i) neg[i] = -data[i];
    return find_peaks_prominence(neg, distance, prominence);
}

// ---------- Zero crossing detection & filtering close ones ----------

std::vector<double> filter_close_zero_crossings(const std::vector<double>& zcrossings, const std::vector<double>& acc, int window=7, int min_separation=10) {
    std::vector<double> filtered;
    if (zcrossings.empty()) return filtered;
    std::vector<double> sorted_z = zcrossings;
    std::sort(sorted_z.begin(), sorted_z.end());
    size_t i = 0;
    while (i < sorted_z.size()) {
        double current = sorted_z[i];
        std::vector<double> cluster = {current};
        size_t j = i + 1;
        while (j < sorted_z.size() && std::abs(sorted_z[j] - current) < min_separation) {
            cluster.push_back(sorted_z[j]);
            ++j;
        }
        if (cluster.size() == 1) {
            filtered.push_back(current);
        } else {
            // choose cluster member with largest local amplitude in window
            std::vector<double> amps;
            for (double z : cluster) {
                int zero_idx = (int)std::round(z);
                int lo = std::max(0, zero_idx - window);
                int hi = std::min((int)acc.size(), zero_idx + window + 1);
                double local_max = -1e300, local_min = 1e300;
                for (int k = lo; k < hi; ++k) { local_max = std::max(local_max, acc[k]); local_min = std::min(local_min, acc[k]); }
                amps.push_back(local_max - local_min);
            }
            int max_idx = 0;
            for (size_t k = 1; k < amps.size(); ++k) if (amps[k] > amps[max_idx]) max_idx = k;
            filtered.push_back(cluster[max_idx]);
        }
        i += cluster.size();
    }
    return filtered;
}

std::vector<double> find_all_zero_crossings(const std::vector<double>& acc, double min_amplitude=0.5, int window=7, double eps=0.05) {
    std::vector<double> zero_crossings;
    int n = acc.size();
    int i = 0;
    while (i < n - 1) {
        if ((acc[i] > eps && acc[i+1] < -eps) || (acc[i] < -eps && acc[i+1] > eps)) {
            double denom = acc[i+1] - acc[i];
            if (std::abs(denom) < 1e-12) { ++i; continue; }
            double z = i - acc[i] / denom;
            int zero_idx = (int)std::round(z);
            int lo = std::max(0, zero_idx - window);
            int hi = std::min(n, zero_idx + window + 1);
            double local_max = -1e300, local_min = 1e300;
            for (int k = lo; k < hi; ++k) { local_max = std::max(local_max, acc[k]); local_min = std::min(local_min, acc[k]); }
            double local_amp = local_max - local_min;
            double slope = std::abs(acc[i+1] - acc[i]);
            if (local_amp >= min_amplitude && slope > 0.15) zero_crossings.push_back(z);
            i += 1;
        } else if (std::abs(acc[i]) <= eps) {
            int start = i;
            while (i < n && std::abs(acc[i]) <= eps) ++i;
            int end = i - 1;
            double mid = 0.5 * (start + end);
            int zero_idx = (int)std::round(mid);
            int lo = std::max(0, zero_idx - window);
            int hi = std::min(n, zero_idx + window + 1);
            double local_max = -1e300, local_min = 1e300;
            for (int k = lo; k < hi; ++k) { local_max = std::max(local_max, acc[k]); local_min = std::min(local_min, acc[k]); }
            double local_amp = local_max - local_min;
            if (local_amp >= min_amplitude) zero_crossings.push_back(mid);
        } else {
            ++i;
        }
    }
    return filter_close_zero_crossings(zero_crossings, acc, window);
}

// ---------- Main processing ----------

struct StrokeCycle {
    int StrokeID;
    int StartZero, PosPeak, MidZero, NegPeak, EndZero;
    double TotalTime, TimeToPosPeak, TimeToNegPeak;
    double Area_Positive, Area_Negative;
};

int main() {
    try {
        std::string RESULT_DIR = "results_full_strokes";
        std::string STROKE_PLOT_DIR = RESULT_DIR + "/stroke_plots";
        fs::create_directories(STROKE_PLOT_DIR);

        // --- Load Data from Excel (sheet name: "輕艇xsens") ---
        xlnt::workbook wb;
        wb.load("data/Canoe xsens dot.xlsx");
        xlnt::worksheet ws;
        try {
            ws = wb.sheet_by_title("輕艇xsens");
        } catch(...) {
            ws = wb.active_sheet();
            std::cout << "Warning: sheet '輕艇xsens' not found; using active sheet.\n";
        }

        // Read columns into vectors by looking up header names
        std::vector<std::string> headers;
        auto row1 = ws.rows(false).begin();
        for (auto cell : *row1) headers.push_back(cell.to_string());
        // find indices of required columns
        auto find_col = [&](const std::string &name)->int{
            for (size_t i=0;i<headers.size();++i) if (headers[i]==name) return i; return -1; };

        int idx_AccX = find_col("Acc_X");
        int idx_AccY = find_col("Acc_Y");
        int idx_GyrX = find_col("Gyr_X");
        int idx_GyrY = find_col("Gyr_Y");
        int idx_GyrZ = find_col("Gyr_Z");
        if (idx_AccY < 0) { std::cerr << "Acc_Y column not found in sheet headers.\n"; return 1; }

        std::vector<double> acc_x, acc_y, gyr_x, gyr_y, gyr_z;
        for (auto it = std::next(ws.rows(false).begin()); it != ws.rows(false).end(); ++it) {
            auto row = *it; 
            // some cells could be missing; guard
            auto cell = row[idx_AccY]; if (!cell.has_value()) break;
            // reading as double (if invalid, treat as 0)
            auto read_val = [&](int col)->double{
                try {
                    return std::stod(row[col].to_string());
                } catch(...) { return 0.0; }
            };
            acc_x.push_back(idx_AccX>=0 ? read_val(idx_AccX) : 0.0);
            acc_y.push_back(read_val(idx_AccY));
            gyr_x.push_back(idx_GyrX>=0 ? read_val(idx_GyrX) : 0.0);
            gyr_y.push_back(idx_GyrY>=0 ? read_val(idx_GyrY) : 0.0);
            gyr_z.push_back(idx_GyrZ>=0 ? read_val(idx_GyrZ) : 0.0);
        }

        // --- Filtering ---
        double fs_rate = 120.0;
        auto acc_y_smooth = lowpass_filter(acc_y, 5.0, fs_rate, 4);

        // --- Peaks ---
        auto pos_peaks = find_peaks_prominence(acc_y_smooth, 5, 0.8);
        auto neg_peaks = find_peaks_prominence_neg(acc_y_smooth, 5, 0.8);

        // --- Zero crossings ---
        auto all_zeros = find_all_zero_crossings(acc_y_smooth, 0.65, 7, 0.1);
        // filter by local amplitude > 0.65
        std::vector<double> all_zeros_filtered;
        for (double z : all_zeros) {
            int zi = (int)std::round(z);
            int lo = std::max(0, zi - 7);
            int hi = std::min((int)acc_y_smooth.size(), zi + 8);
            double local_max = -1e300, local_min = 1e300;
            for (int k = lo; k < hi; ++k) { local_max = std::max(local_max, acc_y_smooth[k]); local_min = std::min(local_min, acc_y_smooth[k]); }
            if (local_max - local_min > 0.65) all_zeros_filtered.push_back(z);
        }

        // --- Full stroke cycles ---
        std::vector<StrokeCycle> stroke_cycles;
        for (size_t i = 0; i + 3 < all_zeros_filtered.size(); ++i) {
            int z1 = (int)std::round(all_zeros_filtered[i]);
            int z2 = (int)std::round(all_zeros_filtered[i+1]);
            int z3 = (int)std::round(all_zeros_filtered[i+2]);
            int z4 = (int)std::round(all_zeros_filtered[i+3]);
            // find pos peak between z1 and z2
            std::vector<int> pos_in_seg;
            for (int p : pos_peaks) if (p > z1 && p < z2) pos_in_seg.push_back(p);
            std::vector<int> neg_in_seg;
            for (int p : neg_peaks) if (p > z2 && p < z3) neg_in_seg.push_back(p);
            if (pos_in_seg.empty() || neg_in_seg.empty()) continue;
            int pos_idx = pos_in_seg[0];
            double maxv = acc_y_smooth[pos_idx];
            for (int p : pos_in_seg) if (acc_y_smooth[p] > maxv) { maxv = acc_y_smooth[p]; pos_idx = p; }
            int neg_idx = neg_in_seg[0];
            double minv = acc_y_smooth[neg_idx];
            for (int p : neg_in_seg) if (acc_y_smooth[p] < minv) { minv = acc_y_smooth[p]; neg_idx = p; }

            double t_total = double(z3 - z1) / fs_rate;
            double t_pos = double(pos_idx - z1) / fs_rate;
            double t_neg = double(neg_idx - z2) / fs_rate;
            // area positive: integrate clipped to >=0 between z1..z2
            std::vector<double> segpos;
            for (int k = z1; k <= z2 && k < (int)acc_y_smooth.size(); ++k) segpos.push_back(std::max(0.0, acc_y_smooth[k]));
            std::vector<double> segneg;
            for (int k = z2; k <= z3 && k < (int)acc_y_smooth.size(); ++k) segneg.push_back(std::min(0.0, acc_y_smooth[k]));
            double pos_area = trapz(segpos, 1.0/fs_rate);
            double neg_area = trapz(segneg, 1.0/fs_rate);

            StrokeCycle sc;
            sc.StrokeID = (int)stroke_cycles.size() + 1;
            sc.StartZero = z1; sc.PosPeak = pos_idx; sc.MidZero = z2; sc.NegPeak = neg_idx; sc.EndZero = z3;
            sc.TotalTime = t_total; sc.TimeToPosPeak = t_pos; sc.TimeToNegPeak = t_neg;
            sc.Area_Positive = pos_area; sc.Area_Negative = neg_area;
            stroke_cycles.push_back(sc);
        }

        // --- Save stroke cycles CSV ---
        fs::create_directories(RESULT_DIR);
        std::ofstream stroke_out(RESULT_DIR + "/stroke_full_cycles.csv");
        stroke_out << "StrokeID,StartZero,PosPeak,MidZero,NegPeak,EndZero,TotalTime,TimeToPosPeak,TimeToNegPeak,Area_Positive,Area_Negative\n";
        for (auto &sc : stroke_cycles) {
            stroke_out << sc.StrokeID << "," << sc.StartZero << "," << sc.PosPeak << "," << sc.MidZero << "," << sc.NegPeak << "," << sc.EndZero << ","
                       << sc.TotalTime << "," << sc.TimeToPosPeak << "," << sc.TimeToNegPeak << "," << sc.Area_Positive << "," << sc.Area_Negative << "\n";
        }
        stroke_out.close();
        std::cout << "Saved " << stroke_cycles.size() << " full stroke cycles\n";

        // --- Plot overview ---
        plt::figure_size(1600, 600);
        plt::plot(acc_y_smooth);
        // plot pos and neg cycle peaks
        std::vector<int> pos_cycle_peaks, neg_cycle_peaks;
        for (auto &sc : stroke_cycles) { pos_cycle_peaks.push_back(sc.PosPeak); neg_cycle_peaks.push_back(sc.NegPeak); }
        std::vector<double> pos_vals, neg_vals;
        for (int p : pos_cycle_peaks) pos_vals.push_back(acc_y_smooth[p]);
        for (int p : neg_cycle_peaks) neg_vals.push_back(acc_y_smooth[p]);
        if (!pos_cycle_peaks.empty()) plt::scatter(pos_cycle_peaks, pos_vals, 30.0);
        if (!neg_cycle_peaks.empty()) plt::scatter(neg_cycle_peaks, neg_vals, 30.0);

        std::vector<int> zeros_cycle;
        for (auto &sc : stroke_cycles) { zeros_cycle.push_back(sc.StartZero); zeros_cycle.push_back(sc.MidZero); zeros_cycle.push_back(sc.EndZero); }
        std::vector<double> zeros_zeros(zeros_cycle.size(), 0.0);
        if (!zeros_cycle.empty()) plt::scatter(zeros_cycle, zeros_zeros, 20.0);

        // fill regions
        for (auto &sc : stroke_cycles) {
            std::vector<double> x_fill;  // x-values
            std::vector<double> y_fill;  // y-values
            for (int k = sc.StartZero; k <= sc.EndZero && k < (int)acc_y_smooth.size(); ++k) {
                x_fill.push_back(static_cast<double>(k));
                y_fill.push_back(acc_y_smooth[k]);
            }
            if (!x_fill.empty()) {
                std::vector<double> y_base(x_fill.size(), 0.0);
                plt::fill_between(x_fill, y_fill, y_base, std::map<std::string,std::string>{}); // <-- fix here
            }
        }

        plt::title("Full Stroke Cycles: 0 → +Peak → 0 → -Peak → 0");
        plt::xlabel("Sample Index"); plt::ylabel("Acc_Y"); plt::grid(true);
        plt::save(RESULT_DIR + "/strokes_full_overview.png");
        plt::close();

        // --- Individual stroke plots ---
        for (size_t i = 0; i < stroke_cycles.size(); ++i) {
            auto &sc = stroke_cycles[i];
            std::vector<int> seg_x; std::vector<double> seg_y;
            for (int k = sc.StartZero; k <= sc.EndZero && k < (int)acc_y_smooth.size(); ++k) { seg_x.push_back(k); seg_y.push_back(acc_y_smooth[k]); }
            if (seg_x.empty()) continue;
            plt::figure_size(1000, 400);
            plt::plot(seg_x, seg_y);
            // horizontal 0
            std::vector<double> zeros_line(seg_x.size(), 0.0);
            plt::plot(seg_x, zeros_line);
            plt::scatter(std::vector<int>{sc.PosPeak}, std::vector<double>{acc_y_smooth[sc.PosPeak]}, 50.0);
            plt::scatter(std::vector<int>{sc.NegPeak}, std::vector<double>{acc_y_smooth[sc.NegPeak]}, 50.0);
            plt::scatter(std::vector<int>{sc.StartZero, sc.MidZero, sc.EndZero}, std::vector<double>{0.0,0.0,0.0}, 30.0);
            plt::title("Full Stroke " + std::to_string(i+1) + ": 0→+Peak→0→-Peak→0");
            plt::xlabel("Sample Index"); plt::ylabel("Acc_Y"); plt::grid(true);
            std::string seg_path = STROKE_PLOT_DIR + "/stroke_full_" + std::to_string(i+1) + ".png";
            plt::save(seg_path);
            plt::close();
        }

        std::cout << "Saved " << stroke_cycles.size() << " individual stroke plots in '" << STROKE_PLOT_DIR << "'\n";

        // --- Phase metrics recalculated per stroke ---
        auto acc_x_smooth = lowpass_filter(acc_x, 5.0, fs_rate, 4);
        auto gyr_x_smooth = lowpass_filter(gyr_x, 5.0, fs_rate, 4);
        auto gyr_y_smooth = lowpass_filter(gyr_y, 5.0, fs_rate, 4);
        auto gyr_z_smooth = lowpass_filter(gyr_z, 5.0, fs_rate, 4);

        struct PhaseMetric { std::map<std::string,double> vals; };
        std::vector<PhaseMetric> phase_metrics;

        for (auto &sc : stroke_cycles) {
            int start_idx = sc.StartZero, pos_idx = sc.PosPeak, mid_idx = sc.MidZero, end_idx = sc.EndZero, neg_idx = sc.NegPeak;
            // Acceleration phase
            std::vector<double> acc_phase;
            for (int k = start_idx; k <= pos_idx && k < (int)acc_y_smooth.size(); ++k) acc_phase.push_back(acc_y_smooth[k]);
            double dur_acc = double(pos_idx - start_idx) / fs_rate;
            double acc_peak = acc_phase.empty() ? 0.0 : *std::max_element(acc_phase.begin(), acc_phase.end());
            double time_to_peak = dur_acc;
            double perc_time_to_peak = 100.0; // per original script
            double RAD = (time_to_peak > 0.0 ? acc_peak / time_to_peak : nan_val);
            double vel_change_acc = trapz(acc_phase, 1.0/fs_rate);
            double mean_abs_ax_acc = mean_abs_diff(slice(acc_x_smooth, start_idx, pos_idx+1));
            double mean_abs_gx_acc = mean_abs_diff(slice(gyr_x_smooth, start_idx, pos_idx+1));
            double mean_abs_gy_acc = mean_abs_diff(slice(gyr_y_smooth, start_idx, pos_idx+1));
            double mean_abs_gz_acc = mean_abs_diff(slice(gyr_z_smooth, start_idx, pos_idx+1));

            // Deceleration phase
            std::vector<double> dec_phase;
            for (int k = pos_idx; k <= end_idx && k < (int)acc_y_smooth.size(); ++k) dec_phase.push_back(acc_y_smooth[k]);
            double dur_dec = double(end_idx - pos_idx) / fs_rate;
            double dec_peak = dec_phase.empty() ? 0.0 : *std::min_element(dec_phase.begin(), dec_phase.end());
            double time_to_neg_peak = double(neg_idx - pos_idx) / fs_rate;
            double perc_time_to_neg = (dur_dec > 0.0 ? 100.0 * time_to_neg_peak / dur_dec : nan_val);
            double vel_change_dec = trapz(dec_phase, 1.0/fs_rate);
            double mean_abs_ax_dec = mean_abs_diff(slice(acc_x_smooth, pos_idx, end_idx+1));
            double mean_abs_gx_dec = mean_abs_diff(slice(gyr_x_smooth, pos_idx, end_idx+1));
            double mean_abs_gy_dec = mean_abs_diff(slice(gyr_y_smooth, pos_idx, end_idx+1));
            double mean_abs_gz_dec = mean_abs_diff(slice(gyr_z_smooth, pos_idx, end_idx+1));

            // Stroke phase
            std::vector<double> gyrx_phase;
            for (int k = start_idx; k <= end_idx && k < (int)gyr_x_smooth.size(); ++k) gyrx_phase.push_back(gyr_x_smooth[k]);
            double gyrx_pos_peak = gyrx_phase.empty() ? 0.0 : *std::max_element(gyrx_phase.begin(), gyrx_phase.end());
            double gyrx_neg_peak = gyrx_phase.empty() ? 0.0 : *std::min_element(gyrx_phase.begin(), gyrx_phase.end());
            double mean_abs_ax_stroke = mean_abs_diff(slice(acc_x_smooth, start_idx, end_idx+1));
            double mean_abs_gx_stroke = mean_abs_diff(slice(gyr_x_smooth, start_idx, end_idx+1));
            double mean_abs_gy_stroke = mean_abs_diff(slice(gyr_y_smooth, start_idx, end_idx+1));
            double mean_abs_gz_stroke = mean_abs_diff(slice(gyr_z_smooth, start_idx, end_idx+1));

            PhaseMetric pm;
            pm.vals["StrokeID"] = sc.StrokeID;
            pm.vals["AccPhase_Duration"] = dur_acc;
            pm.vals["AccPhase_PosPeak"] = acc_peak;
            pm.vals["AccPhase_TimeToPosPeak_%"] = perc_time_to_peak;
            pm.vals["AccPhase_RAD"] = RAD;
            pm.vals["AccPhase_VelChange"] = vel_change_acc;
            pm.vals["AccPhase_MeanAbs_AccX"] = mean_abs_ax_acc;
            pm.vals["AccPhase_MeanAbs_GyrX"] = mean_abs_gx_acc;
            pm.vals["AccPhase_MeanAbs_GyrY"] = mean_abs_gy_acc;
            pm.vals["AccPhase_MeanAbs_GyrZ"] = mean_abs_gz_acc;

            pm.vals["DecPhase_Duration"] = dur_dec;
            pm.vals["DecPhase_NegPeak"] = dec_peak;
            pm.vals["DecPhase_TimeToNegPeak_%"] = perc_time_to_neg;
            pm.vals["DecPhase_VelChange"] = vel_change_dec;
            pm.vals["DecPhase_MeanAbs_AccX"] = mean_abs_ax_dec;
            pm.vals["DecPhase_MeanAbs_GyrX"] = mean_abs_gx_dec;
            pm.vals["DecPhase_MeanAbs_GyrY"] = mean_abs_gy_dec;
            pm.vals["DecPhase_MeanAbs_GyrZ"] = mean_abs_gz_dec;

            pm.vals["Stroke_GyrX_PosPeak"] = gyrx_pos_peak;
            pm.vals["Stroke_GyrX_NegPeak"] = gyrx_neg_peak;
            pm.vals["Stroke_MeanAbs_AccX"] = mean_abs_ax_stroke;
            pm.vals["Stroke_MeanAbs_GyrX"] = mean_abs_gx_stroke;
            pm.vals["Stroke_MeanAbs_GyrY"] = mean_abs_gy_stroke;
            pm.vals["Stroke_MeanAbs_GyrZ"] = mean_abs_gz_stroke;

            phase_metrics.push_back(pm);
        }

        // --- Save phase metrics to Excel with 3 sheets using xlnt ---
        xlnt::workbook out_wb;
        xlnt::worksheet acc_sheet = out_wb.create_sheet(); acc_sheet.title("AccelerationPhase");
        xlnt::worksheet dec_sheet = out_wb.create_sheet(); dec_sheet.title("DecelerationPhase");
        xlnt::worksheet stroke_sheet = out_wb.create_sheet(); stroke_sheet.title("StrokePhase");

        // headers for acceleration
        std::vector<std::string> acc_hdr = {"StrokeID","AccPhase_Duration","AccPhase_PosPeak","AccPhase_TimeToPosPeak_%","AccPhase_RAD","AccPhase_VelChange","AccPhase_MeanAbs_AccX","AccPhase_MeanAbs_GyrX","AccPhase_MeanAbs_GyrY","AccPhase_MeanAbs_GyrZ"};
        std::vector<std::string> dec_hdr = {"StrokeID","DecPhase_Duration","DecPhase_NegPeak","DecPhase_TimeToNegPeak_%","DecPhase_VelChange","DecPhase_MeanAbs_AccX","DecPhase_MeanAbs_GyrX","DecPhase_MeanAbs_GyrY","DecPhase_MeanAbs_GyrZ"};
        std::vector<std::string> stroke_hdr = {"StrokeID","Stroke_GyrX_PosPeak","Stroke_GyrX_NegPeak","Stroke_MeanAbs_AccX","Stroke_MeanAbs_GyrX","Stroke_MeanAbs_GyrY","Stroke_MeanAbs_GyrZ"};

        // write headers
        for (size_t c = 0; c < acc_hdr.size(); ++c) acc_sheet.cell(1, c+1).value(acc_hdr[c]);
        for (size_t c = 0; c < dec_hdr.size(); ++c) dec_sheet.cell(1, c+1).value(dec_hdr[c]);
        for (size_t c = 0; c < stroke_hdr.size(); ++c) stroke_sheet.cell(1, c+1).value(stroke_hdr[c]);

        // write rows
        for (size_t r = 0; r < phase_metrics.size(); ++r) {
            auto &pm = phase_metrics[r].vals;
            int row_no = (int)r + 2;
            // acc sheet
            for (size_t c = 0; c < acc_hdr.size(); ++c) acc_sheet.cell(row_no, c+1).value(pm[acc_hdr[c]]);
            for (size_t c = 0; c < dec_hdr.size(); ++c) dec_sheet.cell(row_no, c+1).value(pm[dec_hdr[c]]);
            for (size_t c = 0; c < stroke_hdr.size(); ++c) stroke_sheet.cell(row_no, c+1).value(pm[stroke_hdr[c]]);
        }

        std::string phase_excel = RESULT_DIR + "/stroke_phase_metrics.xlsx";
        out_wb.save(phase_excel);
        std::cout << "Phase metrics saved to '" << phase_excel << "' with 3 sheets\n";

    } catch (const std::exception &ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 2;
    }

    return 0;
}