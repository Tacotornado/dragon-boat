# ==============================
# Container-Ready Portable Makefile
# Adapted from the known-working build for the engine/src/ layout.
#   Run location (Docker): CWD = engine/  ->  `make all`
#   Sources live under src/, vendored libs under src/libs/.
#   Binary is emitted at engine/stroke_analysis (picked up by the Dockerfile).
# ==============================

CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -static-libstdc++

SRC_DIR := src

# --- Include Directories (vendored project libs under src/libs) ---
# Both libs/cxxopts and libs/cxxopts/include are listed so the header resolves
# whether or not the vendored copy uses an include/ subfolder. Missing -I dirs
# are ignored by g++, so listing both is safe.
INCLUDES := \
  -I$(SRC_DIR) \
  -I$(SRC_DIR)/utilities \
  -I$(SRC_DIR)/libs/eigen \
  -I$(SRC_DIR)/libs/matplotlib-cpp \
  -I$(SRC_DIR)/libs/cxxopts \
  -I$(SRC_DIR)/libs/cxxopts/include \
  -I$(SRC_DIR)/libs/pocketfft \
  $(shell python3 -c "import numpy as c; print('-I' + c.get_include())")

# --- Python Configuration (Dynamic, embed interpreter for matplotlib-cpp) ---
PYTHON_LDFLAGS := $(shell python3-config --embed --ldflags 2>/dev/null || python3-config --ldflags)
PYTHON_CFLAGS := $(shell python3-config --includes)
CXXFLAGS += $(PYTHON_CFLAGS)

# --- Libraries and Linking (xlnt + Python from system paths) ---
LDFLAGS := -lxlnt $(PYTHON_LDFLAGS)

# --- Source Files and Target ---
SRC := \
  $(SRC_DIR)/stroke_analysis.cpp \
  $(SRC_DIR)/utilities/filters.cpp \
  $(SRC_DIR)/utilities/FrequencyAnalyzer.cpp
OBJ := $(SRC:.cpp=.o)
TARGET := stroke_analysis

# ==============================
# Rules
# ==============================

all: $(TARGET)

$(TARGET): $(OBJ)
	@echo "Linking $(TARGET)..."
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@ $(LDFLAGS)

%.o: %.cpp
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	@echo "Cleaning build files..."
	rm -f $(OBJ) $(TARGET)

run: $(TARGET)
	@echo "Running stroke_analysis..."
	./$(TARGET) --input data/Canoe_xsens_dot.csv --output results

debug: CXXFLAGS += -g -O0
debug: clean all

.PHONY: all clean run debug
