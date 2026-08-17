#pragma once

#include <vector>
#include <complex>

namespace Canoe {

// Biquad structure for Direct Form I
struct Biquad {
    double b0 = 1.0, b1 = 0.0, b2 = 0.0;
    double a1 = 0.0, a2 = 0.0;
    
    // State registers
    double x1 = 0.0, x2 = 0.0;
    double y1 = 0.0, y2 = 0.0;

    double process(double x);
};

class CausalAdaptiveFilter {
public:
    CausalAdaptiveFilter(int order = 4, double fs = 120.0, double initial_cutoff = 10.0);

    void update_cutoff(double cutoff);
    double process_sample(double x);
    double get_cutoff() const { return cutoff_; }

private:
    int order_;
    double fs_;
    double cutoff_;
    std::vector<Biquad> biquads_;

    void design_filter(double cutoff);
};

class PsdEstimator {
public:
    PsdEstimator(int window_size = 360, int fft_size = 256, double fs = 120.0);

    // Appends raw sample and returns adaptive cutoff frequency if window is ready and on update step
    bool process_sample(double raw_y, int sample_idx, double& computed_cutoff);

private:
    int window_size_;
    int fft_size_;
    double fs_;
    std::vector<double> buffer_;

    std::vector<double> hann_window_;

    double compute_dominant_frequency();
    void run_fft(std::vector<std::complex<double>>& data);
};

} // namespace Canoe
