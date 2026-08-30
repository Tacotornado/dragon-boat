#include "filters.h"
#include <cmath>

// ==========================================================
// Butterworth Low-pass
// ==========================================================
ButterworthLowpass::ButterworthLowpass(double cutoff_freq, double fs, int order)
    : order_(order), z1_(0.0), z2_(0.0) {
    calculateCoefficients(cutoff_freq, fs);
}

void ButterworthLowpass::calculateCoefficients(double cutoff_freq, double fs) {
    double wc = tan(M_PI * cutoff_freq / fs);
    double k1 = sqrt(2.0) * wc;
    double k2 = wc * wc;
    double a0_inv = 1.0 / (1.0 + k1 + k2);

    a0_ = k2 * a0_inv;
    a1_ = 2.0 * a0_;
    a2_ = a0_;
    b1_ = 2.0 * (k2 - 1.0) * a0_inv;
    b2_ = (1.0 - k1 + k2) * a0_inv;
}

double ButterworthLowpass::process(double x) {
    double y = a0_ * x + a1_ * z1_ + a2_ * z2_ - b1_ * z1_ - b2_ * z2_;
    z2_ = z1_;
    z1_ = y;
    return y;
}

std::vector<double> ButterworthLowpass::apply(const std::vector<double>& data) {
    std::vector<double> output(data.size());
    for (size_t i = 0; i < data.size(); ++i)
        output[i] = process(data[i]);
    return output;
}

// ==========================================================
// Butterworth High-pass
// ==========================================================
ButterworthHighpass::ButterworthHighpass(double cutoff_freq, double fs, int order)
    : order_(order), z1_(0.0), z2_(0.0) {
    calculateCoefficients(cutoff_freq, fs);
}

void ButterworthHighpass::calculateCoefficients(double cutoff_freq, double fs) {
    double wc = tan(M_PI * cutoff_freq / fs);
    double k1 = sqrt(2.0) * wc;
    double k2 = wc * wc;
    double a0_inv = 1.0 / (1.0 + k1 + k2);

    a0_ = 1.0 * a0_inv;
    a1_ = -2.0 * a0_;
    a2_ = a0_;
    b1_ = 2.0 * (k2 - 1.0) * a0_inv;
    b2_ = (1.0 - k1 + k2) * a0_inv;
}

double ButterworthHighpass::process(double x) {
    double y = a0_ * x + a1_ * z1_ + a2_ * z2_ - b1_ * z1_ - b2_ * z2_;
    z2_ = z1_;
    z1_ = y;
    return y;
}

std::vector<double> ButterworthHighpass::apply(const std::vector<double>& data) {
    std::vector<double> output(data.size());
    for (size_t i = 0; i < data.size(); ++i)
        output[i] = process(data[i]);
    return output;
}

// ==========================================================
// Moving Average
// ==========================================================
MovingAverage::MovingAverage(int window_size)
    : window_size_(window_size), buffer_(window_size, 0.0), index_(0), sum_(0.0) {}

double MovingAverage::process(double x) {
    sum_ -= buffer_[index_];
    buffer_[index_] = x;
    sum_ += x;
    index_ = (index_ + 1) % window_size_;
    return sum_ / window_size_;
}

std::vector<double> MovingAverage::apply(const std::vector<double>& data) {
    std::vector<double> output(data.size());
    for (size_t i = 0; i < data.size(); ++i)
        output[i] = process(data[i]);
    return output;
}
