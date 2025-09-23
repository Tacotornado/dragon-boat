// usage sample for filters.h //

#include <iostream>
#include "utilities/filters.h"

int main() {
    // Example data
    std::vector<double> acc_y = {0.1, 0.5, 0.9, 0.7, -0.2, -0.6, -0.9, -0.4, 0.0, 0.3};

    double fs = 120.0;    // Sampling frequency
    double cutoff = 5.0;  // Cutoff frequency in Hz

    // Create filter object
    ButterworthLowpass lp_filter(cutoff, fs);

    // Apply filter
    std::vector<double> acc_y_smooth = lp_filter.apply(acc_y);

    // Print filtered data
    for (double val : acc_y_smooth)
        std::cout << val << std::endl;

    return 0;
}