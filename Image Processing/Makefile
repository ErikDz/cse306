CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2

SRCS = main.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = color_matching

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
clean:
	rm -f $(OBJS) $(TARGET)
