#ifndef FREQUENCY_ANALYZER_H
#define FREQUENCY_ANALYZER_H

#include <vector>

class FrequencyAnalyzer {
public:
    // Only the declaration goes here
    static double estimateStrokeFrequency(const std::vector<double>& data, double fs);
};

#endif