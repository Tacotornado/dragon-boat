# ==============================
# Makefile for stroke_analysis project
# ==============================

# --- Compiler and flags ---
CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -static-libstdc++
# (You can remove -static-libstdc++ if linking statically causes issues)

# --- Include directories ---
INCLUDES := \
  -I/mnt/MPLL/dragonboat/NAS/xlnt-src/include \
  -I/mnt/MPLL/dragonboat/NAS/eigen \
  -I/mnt/MPLL/dragonboat/NAS/matplotlib-cpp \
  -I/mnt/MPLL/dragonboat/miniconda3/include/python3.13 \
  -I/mnt/MPLL/dragonboat/miniconda3/lib/python3.13/site-packages/numpy/core/include\
  -I/mnt/MPLL/dragonboat/dragon-boat/utilities \
  -I/mnt/MPLL/dragonboat/NAS/cxxopts/include

# --- Library directories and linking flags ---
LDFLAGS := \
  -L/mnt/MPLL/dragonboat/NAS/xlnt-src/build/source -lxlnt \
  -L/mnt/MPLL/dragonboat/miniconda3/lib \
  -Wl,-rpath=/mnt/MPLL/dragonboat/NAS/xlnt-src/build/source:/mnt/MPLL/dragonboat/miniconda3/lib \
  $(shell /mnt/MPLL/dragonboat/miniconda3/bin/python3.13-config --embed --ldflags)

# --- Source files and target ---
SRC := stroke_analysis_com.cpp utilities/filters.cpp utilities/auto_cutoff.cpp utilities/fcopt.cpp
OBJ := $(SRC:.cpp=.o)
TARGET := stroke_analysis

# ==============================
# Rules
# ==============================

all: $(TARGET)

# Link the target
$(TARGET): $(OBJ)
	@echo "Linking $(TARGET)..."
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@ $(LDFLAGS)

# Compile each .cpp to .o
%.o: %.cpp
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean build artifacts
clean:
	@echo "Cleaning build files..."
	rm -f $(OBJ) $(TARGET)

# Run program (example usage)
run: $(TARGET)
	@echo "Running stroke_analysis..."
	./$(TARGET) --input data/Canoe_xsens_dot.csv --output results_com

# Debug build (no optimization, with symbols)
debug: CXXFLAGS += -g -O0
debug: clean all

.PHONY: all clean run debug
