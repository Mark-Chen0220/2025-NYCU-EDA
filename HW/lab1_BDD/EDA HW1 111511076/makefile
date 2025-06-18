# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -std=c++14

# Target
TARGET = Lab1

# Source files
SRC = BDD.cpp

# Default target
all: $(TARGET)

# Build rule (no .o files)
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Clean target
clean:
	rm -f $(TARGET)
