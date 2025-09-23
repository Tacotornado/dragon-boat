CXX = g++
CXXFLAGS = -Wall -O2 -Iutilities

# Source files
SRC = filter_sample.cpp utilities/filters.cpp
OBJ = $(SRC:.cpp=.o)

# Target
TARGET = filter_sample

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
