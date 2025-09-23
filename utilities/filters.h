#ifndef FILTERS_H
#define FILTERS_H

#include <vector>

class ButterworthLowpass {
public:
    ButterworthLowpass(double cutoff_freq, double fs, int order = 2);

    // Apply filter to a single sample
    double process(double x);

    // Apply filter to a vector of samples
    std::vector<double> apply(const std::vector<double>& data);

private:
    int order_;
    double a0_, a1_, a2_, b1_, b2_;
    double z1_, z2_;
    void calculateCoefficients(double cutoff_freq, double fs);
};

#endif // FILTERS_H