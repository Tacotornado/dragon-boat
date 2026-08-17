# ==============================
# Makefile for Canoe IMU Analytics C++ Computational Engine
# ==============================

CXX := g++
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra

# Include directories
INCLUDES := -Isrc

# Source files
SRC := src/Filter.cpp src/StrokeSegmenter.cpp src/main.cpp
OBJ := $(SRC:.cpp=.o)
TARGET := stroke_analyzer_cpp

all: $(TARGET)

$(TARGET): $(OBJ)
	@echo "Linking $(TARGET)..."
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@

%.o: %.cpp
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	@echo "Cleaning build files..."
	rm -f $(OBJ) $(TARGET)

run: $(TARGET)
	@echo "Running $(TARGET)..."
	./$(TARGET)

.PHONY: all clean run
