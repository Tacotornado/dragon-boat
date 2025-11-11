// signal_features.h
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

// ---- Peak detection (positive, refined with interpolation) ----
inline std::vector<double> find_peaks_prominence(const std::vector<double> &sig, int distance, double prominence) {
    std::vector<double> peaks;
    int n = (int)sig.size();

    for (int i = distance; i < n - distance; i++) {
        if (sig[i] > sig[i - 1] && sig[i] > sig[i + 1]) {
            bool ok = true;
            for (int j = 1; j <= distance; j++) {
                if (sig[i] < sig[i - j] || sig[i] < sig[i + j]) { ok = false; break; }
            }
            if (!ok) continue;

            double left_min = sig[i];
            for (int j = i; j >= 0; j--) left_min = std::min(left_min, sig[j]);
            double right_min = sig[i];
            for (int j = i; j < n; j++) right_min = std::min(right_min, sig[j]);
            double prom = sig[i] - std::max(left_min, right_min);
            if (prom < prominence) continue;

            // --- Parabolic interpolation for sub-sample precision ---
            double denom = (sig[i - 1] - 2 * sig[i] + sig[i + 1]);
            double frac_offset = 0.0;
            if (std::abs(denom) > 1e-12)
                frac_offset = 0.5 * (sig[i - 1] - sig[i + 1]) / denom;
            peaks.push_back((double)i + frac_offset);
        }
    }
    return peaks;
}

// ---- Peak detection (negative, refined with interpolation) ----
inline std::vector<double> find_peaks_prominence_neg(const std::vector<double> &sig, int distance, double prominence) {
    std::vector<double> peaks;
    int n = (int)sig.size();

    for (int i = distance; i < n - distance; i++) {
        if (sig[i] < sig[i - 1] && sig[i] < sig[i + 1]) {
            bool ok = true;
            for (int j = 1; j <= distance; j++) {
                if (sig[i] > sig[i - j] || sig[i] > sig[i + j]) { ok = false; break; }
            }
            if (!ok) continue;

            double left_max = sig[i];
            for (int j = i; j >= 0; j--) left_max = std::max(left_max, sig[j]);
            double right_max = sig[i];
            for (int j = i; j < n; j++) right_max = std::max(right_max, sig[j]);
            double prom = std::min(left_max, right_max) - sig[i];
            if (prom < prominence) continue;

            // --- Parabolic interpolation for sub-sample precision ---
            double denom = (sig[i - 1] - 2 * sig[i] + sig[i + 1]);
            double frac_offset = 0.0;
            if (std::abs(denom) > 1e-12)
                frac_offset = 0.5 * (sig[i - 1] - sig[i + 1]) / denom;
            peaks.push_back((double)i + frac_offset);
        }
    }
    return peaks;
}

// ---- Linear interpolation helper by fractional index ----
inline double interp_at(const std::vector<double>& sig, double pos) {
    if (sig.empty()) return 0.0;
    int n = (int)sig.size();
    if (pos <= 0.0) return sig.front();
    if (pos >= n - 1) return sig.back();
    int i0 = (int)std::floor(pos);
    int i1 = i0 + 1;
    double frac = pos - (double)i0;
    return sig[i0] * (1.0 - frac) + sig[i1] * frac;
}

// ---- Zero crossing detection (returns fractional indices; MATLAB-style) ----
inline std::vector<double> find_all_zero_crossings(
    const std::vector<double> &sig,
    double min_amp,    // minimum local amplitude (peak-to-trough) in the window
    int window,        // half-window for local amplitude check (samples)
    double slope_thr,  // minimum slope magnitude between adjacent samples
    int min_spacing,   // minimum spacing (samples) between consecutive zeros
    double hysteresis) // kept for API parity (unused here)
{
    std::vector<double> zeros;
    int n = (int)sig.size();
    int last_z = -min_spacing;

    for (int i = 1; i < n; ++i) {
        double y1 = sig[i - 1];
        double y2 = sig[i];

        // require sign change between adjacent samples
        bool sign_change = (y1 <= 0.0 && y2 > 0.0) || (y1 >= 0.0 && y2 < 0.0);
        if (!sign_change) continue;

        // local amplitude check
        int lo = std::max(0, i - window);
        int hi = std::min(n - 1, i + window);
        double local_max = *std::max_element(sig.begin() + lo, sig.begin() + hi + 1);
        double local_min = *std::min_element(sig.begin() + lo, sig.begin() + hi + 1);
        if ((local_max - local_min) < min_amp) continue;

        // slope check
        double slope = y2 - y1;
        if (std::abs(slope) < slope_thr) continue;

        // spacing check (relative to last integer crossing index)
        if (i - last_z < min_spacing) continue;
        last_z = i;

        // linear interpolation to compute fractional index of zero crossing
        double denom = (y2 - y1);
        if (std::abs(denom) > 1e-12) {
            double frac = -y1 / denom;                // fraction along [i-1, i]
            double zero_pos = (double)(i - 1) + frac; // fractional sample index
            zeros.push_back(zero_pos);
        } else {
            // fallback to integer index if denom is zero (extremely rare)
            zeros.push_back((double)i);
        }
    }
    return zeros;
}
