//math_utils.h
#pragma once
#include <vector>
#include <numeric>

inline double trapz(const std::vector<double> &y, double dx) {
    if (y.size() < 2) return 0.0;
    double s = 0.0;
    for (size_t i = 1; i < y.size(); i++)
        s += 0.5 * (y[i-1] + y[i]) * dx;
    return s;
}

// ---- MATLAB-like zero crossing (returns fractional x positions) ----
inline std::vector<double> find_zero_crossings_interp(const std::vector<double>& x,
                                                     const std::vector<double>& y) {
    std::vector<double> zeros;
    if (x.size() != y.size() || x.size() < 2)
        return zeros;

    for (size_t i = 1; i < x.size(); ++i) {
        double y1 = y[i - 1];
        double y2 = y[i];
        if ((y1 > 0 && y2 < 0) || (y1 < 0 && y2 > 0)) {
            double frac = -y1 / (y2 - y1);
            double x0 = x[i - 1] + frac * (x[i] - x[i - 1]);
            zeros.push_back(x0);
        }
    }
    return zeros;
}

inline std::vector<double> interp1(const std::vector<double>& x,
                                  const std::vector<double>& y,
                                  const std::vector<double>& xi) {
    std::vector<double> yi;
    if (x.size() != y.size() || x.size() < 2) return yi;
    for (double x0 : xi) {
        if (x0 <= x.front()) { yi.push_back(y.front()); continue; }
        if (x0 >= x.back())  { yi.push_back(y.back());  continue; }

        auto it = std::upper_bound(x.begin(), x.end(), x0);
        size_t i = std::distance(x.begin(), it);
        double x1 = x[i - 1], x2 = x[i];
        double y1 = y[i - 1], y2 = y[i];
        double frac = (x0 - x1) / (x2 - x1);
        yi.push_back(y1 + frac * (y2 - y1));
    }
    return yi;
}



