#ifndef AUTO_CUTOFF_H
#define AUTO_CUTOFF_H

#include <vector>

namespace AutoCutoff {
    // estimate cutoff using spectrum energy knee/elbow method
    // signal: time-domain samples
    // fs: sampling frequency (Hz)
    // returns cutoff frequency in Hz
    double estimateCutoffKneePoint(const std::vector<double>& signal, double fs,
                                   double f_min = 3.0, double f_max = 10.0,
                                   double fallback = 5.0);
}

#endif // AUTO_CUTOFF_H