# Makefile for stroke_analysis project

# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17 -static-libstdc++

# Include directories
INCLUDES := \
  -I/mnt/MPLL/dragonboat/NAS/xlnt-src/include \
  -I/mnt/MPLL/dragonboat/NAS/eigen \
  -I/mnt/MPLL/dragonboat/NAS/matplotlib-cpp \
  -I/mnt/MPLL/dragonboat/miniconda3/include/python3.13 \
  -I/mnt/MPLL/dragonboat/miniconda3/lib/python3.13/site-packages/numpy/_core/include \
  -I/mnt/MPLL/dragonboat/dragon-boat/utilities

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
$(TARGET):
	$(CXX) $(CXXFLAGS) $(INCLUDES) stroke_analysis.cpp utilities/filters.cpp -o $@ $(LDFLAGS)

# Compile source files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(OBJ) $(TARGET)

# Run the program
run: $(TARGET)
	./$(TARGET)
