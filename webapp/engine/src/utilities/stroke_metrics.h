// stroke_metrics.h
#pragma once
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "math_utils.h"

// compute positive & negative areas between two zero crossings
inline std::pair<double,double> computeAreas(const std::vector<double> &sig, int z1, int z2, double fs) {
    std::vector<double> pos, neg;
    pos.reserve(z2 - z1 + 1);
    neg.reserve(z2 - z1 + 1);

    for (int k = z1; k <= z2; k++) {
        if (sig[k] > 0) pos.push_back(sig[k]); else pos.push_back(0.0);
        if (sig[k] < 0) neg.push_back(sig[k]); else neg.push_back(0.0);
    }
    double area_pos = trapz(pos, 1.0/fs);
    double area_neg = std::abs(trapz(neg, 1.0/fs));
    return {area_pos, area_neg};
}

// compute time metrics
struct StrokeTiming {
    double totalTime;
    double timeToPosPeak;
    double timeToNegPeak;
};

inline StrokeTiming computeTiming(int z1, int posPeak, int midZero, int negPeak, int z3, double fs) {
    StrokeTiming t;
    t.totalTime   = (z3 - z1) / fs;
    t.timeToPosPeak = (posPeak - z1) / fs;
    t.timeToNegPeak = (negPeak - midZero) / fs;
    return t;
}
