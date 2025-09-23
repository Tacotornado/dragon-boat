#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <algorithm>
#include <vector>

// Set up //
std::string DATA_DIR = "sensor_data.csv";
std::string RESULT_DIR = "results_full_strokes";

int main()
{
    std::ifstream data_file(DATA_DIR);
    if (!data_file.is_open()) {
        std::cerr << "Error opening data file: " << DATA_DIR << std::endl;
        return 1;
    }
}