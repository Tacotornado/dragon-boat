#ifndef FILTERS_H
#define FILTERS_H

#include <vector>

// -------------------- Base Filter Class --------------------
class IFilter {
public:
    virtual ~IFilter() = default;
    virtual double process(double x) = 0;                     // process single sample
    virtual std::vector<double> apply(const std::vector<double>& data) = 0; // process vector
};

// -------------------- Butterworth Low-pass --------------------
class ButterworthLowpass : public IFilter {
public:
    ButterworthLowpass(double cutoff_freq, double fs, int order = 2);
    double process(double x) override;
    std::vector<double> apply(const std::vector<double>& data) override;

private:
    int order_;
    double a0_, a1_, a2_, b1_, b2_;
    double z1_, z2_;
    void calculateCoefficients(double cutoff_freq, double fs);
};

// -------------------- Butterworth High-pass --------------------
class ButterworthHighpass : public IFilter {
public:
    ButterworthHighpass(double cutoff_freq, double fs, int order = 2);
    double process(double x) override;
    std::vector<double> apply(const std::vector<double>& data) override;

private:
    int order_;
    double a0_, a1_, a2_, b1_, b2_;
    double z1_, z2_;
    void calculateCoefficients(double cutoff_freq, double fs);
};

// -------------------- Moving Average Filter --------------------
class MovingAverage : public IFilter {
public:
    MovingAverage(int window_size);
    double process(double x) override;
    std::vector<double> apply(const std::vector<double>& data) override;

private:
    int window_size_;
    std::vector<double> buffer_;
    int index_;
    double sum_;
};

#endif // FILTERS_H
