#include "filters.h"
#include <cmath>

ButterworthLowpass::ButterworthLowpass(double cutoff_freq, double fs, int order)
    : order_(order), z1_(0.0), z2_(0.0)
{
    calculateCoefficients(cutoff_freq, fs);
}

void ButterworthLowpass::calculateCoefficients(double cutoff_freq, double fs) {
    // Pre-warped cutoff for bilinear transform
    double wc = tan(M_PI * cutoff_freq / fs);
    double k1 = sqrt(2.0) * wc;
    double k2 = wc * wc;
    double a0_inv = 1.0 / (1.0 + k1 + k2);

    a0_ = k2 * a0_inv;
    a1_ = 2.0 * a0_;
    a2_ = a0_;
    b1_ = 2.0 * (k2 - 1.0) * a0_inv;
    b2_ = (1.0 - k1 + k2) * a0_inv;
    z1_ = z2_ = 0.0;
}

double ButterworthLowpass::process(double x) {
    double y = a0_ * x + a1_ * z1_ + a2_ * z2_ - b1_ * z1_ - b2_ * z2_;
    z2_ = z1_;
    z1_ = y;
    return y;
}

std::vector<double> ButterworthLowpass::apply(const std::vector<double>& data) {
    std::vector<double> output(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        output[i] = process(data[i]);
    }
    return output;
}