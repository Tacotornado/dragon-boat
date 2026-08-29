# ==============================
# Container-Ready Portable Makefile
# ==============================

CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -static-libstdc++

# --- Include Directories (Using local project libs) ---
INCLUDES := \
  -Iutilities \
  -Ilibs/eigen \
  -Ilibs/matplotlib-cpp \
  -Ilibs/cxxopts/include \
  -Ilibs/pocketfft \
  $(shell python3 -c "import numpy as c; print('-I' + c.get_include())")

# --- Python Configuration (Dynamic) ---
PYTHON_LDFLAGS := $(shell python3-config --embed --ldflags 2>/dev/null || python3-config --ldflags)
PYTHON_CFLAGS := $(shell python3-config --includes)
CXXFLAGS += $(PYTHON_CFLAGS)

# --- Libraries and Linking ---
# xlnt and Python are linked via standard system paths inside the container
LDFLAGS := -lxlnt $(PYTHON_LDFLAGS)

# --- Source Files and Target ---
SRC := stroke_analysis.cpp utilities/filters.cpp utilities/FrequencyAnalyzer.cpp
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