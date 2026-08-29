#include "FrequencyAnalyzer.h"
#include "pocketfft_hdronly.h"
#include <complex>
#include <vector>
#include <algorithm>

double FrequencyAnalyzer::estimateStrokeFrequency(const std::vector<double>& data, double fs) {
    size_t n = data.size();
    std::vector<std::complex<double>> buffer(n);

    for (size_t i = 0; i < n; ++i) buffer[i] = {data[i], 0.0};

    pocketfft::shape_t shape{n};
    pocketfft::stride_t stride{sizeof(std::complex<double>)};
    pocketfft::shape_t axes{0};
    pocketfft::c2c(shape, stride, stride, axes, true, buffer.data(), buffer.data(), 1.0);

    double max_val = -1.0;
    int peak_idx = 0;
    for(size_t i = 1; i < n/2; ++i) {
        double psd_val = std::norm(buffer[i]);
        if(psd_val > max_val) {
            max_val = psd_val;
            peak_idx = i;
        }
    }
    return (peak_idx * fs) / n;
}

/*
Usage example

1. Get dominant cadence from last 512 samples
double f_stroke = FrequencyAnalyzer::estimateStrokeFrequency(sensor.acc_z, 120.0);

2. Set filter cutoff to 2x the stroke frequency
double f_cutoff = f_stroke * 2.0;

3. Initialize existing Butterworth filter
ButterworthLowpass filter(f_cutoff, 120.0, 2);

*/