CXX        = g++
CXXFLAGS   = -O3 -Ofast -std=c++20 -fopenmp
LDFLAGS    = -fopenmp

SOURCES    = $(wildcard *.cpp)
OBJECTS    = $(SOURCES:.cpp=.o)
TARGET     = joint_fibre_optimizer

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean
