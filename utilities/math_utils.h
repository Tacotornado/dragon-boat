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
