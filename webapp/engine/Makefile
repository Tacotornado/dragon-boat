# Builds the stroke_analysis binary that the Next.js server spawns.
#
# Expected layout (copy these in from your C++ repo):
#   engine/src/stroke_analysis.cpp   <- main program (already here)
#   engine/utilities/*.cpp *.h *.hpp <- filters.cpp, metrics.cpp, filters.h,
#                                       signal_features.h, stroke_metrics.h,
#                                       math_utils.h, cxxopts.hpp
#   engine/vendor/matplotlibcpp.h    <- header-only matplotlib-cpp
#
# Build:  make -C engine
# Output: engine/stroke_analysis
# Run:    engine/stroke_analysis -i input.csv -o results_dir

CXX      ?= g++
TARGET   := stroke_analysis
SRC      := src/stroke_analysis.cpp $(wildcard utilities/*.cpp)

PYTHON_INCLUDES := $(shell python3-config --includes)
PYTHON_LIBS     := $(shell python3-config --ldflags --embed 2>/dev/null || python3-config --ldflags)
NUMPY_INCLUDE   := $(shell python3 -c "import numpy; print(numpy.get_include())" 2>/dev/null)

CXXFLAGS := -std=c++17 -O2 -Iutilities -Ivendor -Isrc $(PYTHON_INCLUDES)
ifneq ($(NUMPY_INCLUDE),)
CXXFLAGS += -I$(NUMPY_INCLUDE)
endif

LDLIBS := $(PYTHON_LIBS) -lxlnt

.PHONY: all clean
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDLIBS)

clean:
	rm -f $(TARGET)
