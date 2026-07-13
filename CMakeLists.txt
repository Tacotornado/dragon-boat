cmake_minimum_required(VERSION 3.12)
project(CanoeImuAnalytics CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add source files
set(SOURCES
    src/Filter.cpp
    src/StrokeSegmenter.cpp
    src/main.cpp
)

# Header search paths
include_directories(src)

# Define executable
add_executable(stroke_analyzer_cpp ${SOURCES})
