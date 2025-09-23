#include <iostream>
#include "filters.h"

int main() {
    std::vector<double> acc_y = {0.1, 0.5, 0.9, 0.7, -0.2, -0.6, -0.9, -0.4, 0.0, 0.3};
    double fs = 120.0;
    double cutoff = 5.0;

    // --- Low-pass filter ---
    ButterworthLowpass lp(cutoff, fs);
    auto lowpass_result = lp.apply(acc_y);

    std::cout << "Low-pass:\n";
    for (auto v : lowpass_result) std::cout << v << " ";
    std::cout << "\n";

    // --- High-pass filter ---
    ButterworthHighpass hp(cutoff, fs);
    auto highpass_result = hp.apply(acc_y);

    std::cout << "High-pass:\n";
    for (auto v : highpass_result) std::cout << v << " ";
    std::cout << "\n";

    // --- Moving average filter ---
    MovingAverage ma(3);
    auto ma_result = ma.apply(acc_y);

    std::cout << "Moving average:\n";
    for (auto v : ma_result) std::cout << v << " ";
    std::cout << "\n";

    return 0;
}
