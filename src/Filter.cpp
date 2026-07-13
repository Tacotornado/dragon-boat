#include "Filter.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>

namespace Canoe {

const double PI = 3.14159265358979323846;

// --- Biquad Implementation ---
double Biquad::process(double x) {
    double y = b0 * x + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
    x2 = x1;
    x1 = x;
    y2 = y1;
    y1 = y;
    return y;
}

// --- CausalAdaptiveFilter Implementation ---
CausalAdaptiveFilter::CausalAdaptiveFilter(int order, double fs, double initial_cutoff)
    : order_(order), fs_(fs), cutoff_(initial_cutoff) {
    int num_sections = (order + 1) / 2;
    biquads_.resize(num_sections);
    design_filter(initial_cutoff);
}

void CausalAdaptiveFilter::update_cutoff(double cutoff) {
    if (std::abs(cutoff_ - cutoff) > 1e-5) {
        cutoff_ = cutoff;
        design_filter(cutoff);
    }
}

void CausalAdaptiveFilter::design_filter(double cutoff) {
    double theta = PI * cutoff / fs_;
    double tan_val = std::tan(theta);
    if (tan_val == 0.0) return;
    double c = 1.0 / tan_val;
    double c2 = c * c;

    int num_sections = biquads_.size();
    std::vector<double> norms(num_sections);

    for (int s = 0; s < num_sections; ++s) {
        double g = 0.0;
        if (order_ == 2) {
            g = std::sqrt(2.0);
        } else if (order_ == 4) {
            if (s == 0) {
                g = 2.0 * std::cos(1.0 * PI / 8.0); // SciPy Section 0 (damping ~1.84775907)
            } else {
                g = 2.0 * std::cos(3.0 * PI / 8.0); // SciPy Section 1 (damping ~0.76536686)
            }
        } else {
            double angle = PI * (2.0 * (order_ / 2 - 1 - s) + 1.0) / (2.0 * order_);
            g = 2.0 * std::cos(PI / 2.0 - angle);
        }

        double norm = c2 + g * c + 1.0;
        norms[s] = norm;
        biquads_[s].a1 = 2.0 * (1.0 - c2) / norm;
        biquads_[s].a2 = (c2 - g * c + 1.0) / norm;
    }

    // Compute overall gain and assign to first section, keeping other sections normalized
    double g_overall = 1.0;
    for (int s = 0; s < num_sections; ++s) {
        g_overall *= (1.0 / norms[s]);
    }

    for (int s = 0; s < num_sections; ++s) {
        if (s == 0) {
            biquads_[s].b0 = g_overall;
            biquads_[s].b1 = 2.0 * g_overall;
            biquads_[s].b2 = g_overall;
        } else {
            biquads_[s].b0 = 1.0;
            biquads_[s].b1 = 2.0;
            biquads_[s].b2 = 1.0;
        }
    }
}

double CausalAdaptiveFilter::process_sample(double x) {
    double val = x;
    for (auto& bq : biquads_) {
        val = bq.process(val);
    }
    return val;
}

// --- PsdEstimator Implementation ---
PsdEstimator::PsdEstimator(int window_size, int fft_size, double fs)
    : window_size_(window_size), fft_size_(fft_size), fs_(fs) {
    buffer_.reserve(window_size_);
    
    // Precompute Hann window of size fft_size_ (256) to match SciPy's nperseg=256
    hann_window_.resize(fft_size_);
    for (int i = 0; i < fft_size_; ++i) {
        hann_window_[i] = 0.5 * (1.0 - std::cos(2.0 * PI * i / (fft_size_ - 1)));
    }
}

bool PsdEstimator::process_sample(double raw_y, int sample_idx, double& computed_cutoff) {
    if ((int)buffer_.size() < window_size_) {
        buffer_.push_back(raw_y);
    } else {
        for (size_t idx = 0; idx < buffer_.size() - 1; ++idx) {
            buffer_[idx] = buffer_[idx + 1];
        }
        buffer_[window_size_ - 1] = raw_y;
    }

    bool is_full = ((int)buffer_.size() == window_size_);
    bool is_update_step = (sample_idx >= window_size_ && (sample_idx - window_size_) % 120 == 0);

    if (is_full && is_update_step) {
        double f_dom = compute_dominant_frequency();
        double target_cutoff = 5.0 * f_dom;
        computed_cutoff = std::max(5.0, std::min(20.0, target_cutoff));
        return true;
    }

    return false;
}

double PsdEstimator::compute_dominant_frequency() {
    // 1. Extract first 256 samples from the 360 window (indices 0 to 255)
    double sum = 0.0;
    for (int i = 0; i < fft_size_; ++i) {
        sum += buffer_[i];
    }
    double mean = sum / fft_size_;

    // 2. Subtract mean (detrend constant) and apply Hann window
    std::vector<std::complex<double>> fft_data(fft_size_, 0.0);
    for (int i = 0; i < fft_size_; ++i) {
        fft_data[i] = (buffer_[i] - mean) * hann_window_[i];
    }

    // 3. Run 256-point FFT
    run_fft(fft_data);

    // 4. Find dominant frequency peak in [0.3, 2.0] Hz
    double bin_width = fs_ / fft_size_; // 120.0 / 256 = 0.46875 Hz
    
    int k_min = (int)std::ceil(0.3 / bin_width);   // ceil(0.64) = 1
    int k_max = (int)std::floor(2.0 / bin_width);  // floor(4.26) = 4
    
    k_min = std::max(0, k_min);
    k_max = std::min(fft_size_ / 2, k_max);

    int max_k = k_min;
    double max_power = -1.0;

    for (int k = k_min; k <= k_max; ++k) {
        double power = std::norm(fft_data[k]);
        if (power > max_power) {
            max_power = power;
            max_k = k;
        }
    }

    return max_k * bin_width;
}

static unsigned int reverse_bits(unsigned int val, int width) {
    unsigned int result = 0;
    for (int i = 0; i < width; i++) {
        result = (result << 1) | (val & 1);
        val >>= 1;
    }
    return result;
}

void PsdEstimator::run_fft(std::vector<std::complex<double>>& data) {
    int n = data.size();
    int log2n = 0;
    while ((1 << log2n) < n) log2n++;

    for (int i = 0; i < n; i++) {
        int rev = reverse_bits(i, log2n);
        if (i < rev) {
            std::swap(data[i], data[rev]);
        }
    }

    for (int len = 2; len <= n; len <<= 1) {
        double angle = -2.0 * PI / len;
        std::complex<double> wlen(std::cos(angle), std::sin(angle));
        for (int i = 0; i < n; i += len) {
            std::complex<double> w(1.0, 0.0);
            for (int j = 0; j < len / 2; j++) {
                std::complex<double> u = data[i + j];
                std::complex<double> v = data[i + j + len / 2] * w;
                data[i + j] = u + v;
                data[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

} // namespace Canoe
