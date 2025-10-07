#pragma once
#include <bits/stdc++.h>
#include "math_utils.h"

inline double mean(const std::vector<double> &v) {
    if (v.empty()) return 0.0;
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

inline double mean_abs_diff(const std::vector<double> &v) {
    if (v.size() < 2) return 0.0;
    double sum = 0.0;
    for (size_t i = 1; i < v.size(); i++) sum += std::abs(v[i] - v[i-1]);
    return sum / (v.size()-1);
}

// ---- Peak detection (positive) ----
inline std::vector<int> find_peaks_prominence(const std::vector<double> &sig, int distance, double prominence) {
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

// ---- Peak detection (negative) ----
inline std::vector<int> find_peaks_prominence_neg(const std::vector<double> &sig, int distance, double prominence) {
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
inline std::vector<int> find_all_zero_crossings(const std::vector<double> &sig,
                                         double min_amp, int window,
                                         double slope_thr, int min_spacing,
                                         double hysteresis) {
    std::vector<int> zeros;
    int last_z = -min_spacing;
    auto sign = [&](double x){ return (x > hysteresis) - (x < -hysteresis); };
    int n = (int)sig.size();

    for (int i = 1; i < n; i++) {
        int s_prev = sign(sig[i-1]);
        int s_curr = sign(sig[i]);
        if (s_prev == 0) s_prev = sign(sig[ std::max(0, i-2) ]);
        if (s_prev != 0 && s_curr != 0 && s_prev != s_curr) {
            double slope = sig[i] - sig[i-1];
            int lo = std::max(0, i-window);
            int hi = std::min(n-1, i+window);
            double local_max = *std::max_element(sig.begin()+lo, sig.begin()+hi+1);
            double local_min = *std::min_element(sig.begin()+lo, sig.begin()+hi+1);
            if ((local_max - local_min) > min_amp && std::abs(slope) > slope_thr) {
                if (i - last_z >= min_spacing) {
                    // refine: find sample closest to zero
                    int best_idx = lo;
                    double best_abs = std::abs(sig[lo]);
                    for (int k = lo+1; k <= hi; ++k) {
                        double a = std::abs(sig[k]);
                        if (a < best_abs) { best_abs = a; best_idx = k; }
                    }
                    zeros.push_back(best_idx);
                    last_z = i;
                }
            }
        }
    }
    return zeros;
}
