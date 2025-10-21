<<<<<<< HEAD
# ==============================
# Makefile for stroke_analysis project
# ==============================

# --- Compiler and flags ---
CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -static-libstdc++
# (You can remove -static-libstdc++ if linking statically causes issues)

# --- Include directories ---
=======
# Makefile for stroke_analysis project

# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17 -static-libstdc++

# Include directories
>>>>>>> f46fcdd (makefile update)
INCLUDES := \
  -I/mnt/MPLL/dragonboat/NAS/xlnt-src/include \
  -I/mnt/MPLL/dragonboat/NAS/eigen \
  -I/mnt/MPLL/dragonboat/NAS/matplotlib-cpp \
  -I/mnt/MPLL/dragonboat/miniconda3/include/python3.13 \
<<<<<<< HEAD
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
SRC := stroke_analysis_fcopt_d2.cpp utilities/filters.cpp utilities/auto_cutoff.cpp utilities/fcopt.cpp
OBJ := $(SRC:.cpp=.o)
TARGET := stroke_analysis# ==============================
# Rules
# ==============================
=======
  -I/mnt/MPLL/dragonboat/miniconda3/lib/python3.13/site-packages/numpy/_core/include \
  -I/mnt/MPLL/dragonboat/dragon-boat/utilities
>>>>>>> f46fcdd (makefile update)

# Library directories and linking
LDFLAGS := \
  -L/mnt/MPLL/dragonboat/NAS/xlnt-src/build/source -lxlnt \
  -L/mnt/MPLL/dragonboat/miniconda3/lib \
  -Wl,-rpath=/mnt/MPLL/dragonboat/NAS/xlnt-src/build/source:/mnt/MPLL/dragonboat/miniconda3/lib \
  $(shell /mnt/MPLL/dragonboat/miniconda3/bin/python3.13-config --embed --ldflags)

# Source files and output
SRC := stroke_analysis.cpp utilities/filters.cpp
OBJ := $(SRC:.cpp=.o)
TARGET := stroke_analysis

# Default rule
all: $(TARGET)

# Link the target
$(TARGET): $(OBJ)
<<<<<<< HEAD
	# --- Compiler and flags ---
CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -static-libstdc++

# --- Include directories ---
INCLUDES := \
  -I/mnt/MPLL/dragonboat/NAS/xlnt-src/include \
  -I/mnt/MPLL/dragonboat/NAS/eigen \
  -I/mnt/MPLL/dragonboat/NAS/matplotlib-cpp \
  -I/mnt/MPLL/dragonboat/miniconda3/include/python3.13 \
  -I/mnt/MPLL/dragonboat/miniconda3/lib/python3.13/site-packages/numpy/core/include \
  -I/mnt/MPLL/dragonboat/dragon-boat/utilities \
  -I/mnt/MPLL/dragonboat/NAS/cxxopts/include

# --- Library directories and linking flags ---
LDFLAGS := \
  -L/mnt/MPLL/dragonboat/NAS/xlnt-src/build/source -lxlnt \
  -L/mnt/MPLL/dragonboat/miniconda3/lib \
  -Wl,-rpath=/mnt/MPLL/dragonboat/NAS/xlnt-src/build/source:/mnt/MPLL/dragonboat/miniconda3/lib \
  $(shell /mnt/MPLL/dragonboat/miniconda3/bin/python3.13-config --embed --ldflags)

# --- Source files and target ---
# Note: Ensure these files actually exist in your directory
SRC := stroke_analysis_fcopt_d2.cpp utilities/filters.cpp utilities/auto_cutoff.cpp utilities/fcopt.cpp
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
	./$(TARGET)

debug: CXXFLAGS += -g -O0
debug: clean all

.PHONY: all clean run debug
