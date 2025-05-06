CXX = c++
CXXFLAGS = -std=c++17 -O3 -Wall
LDFLAGS =
LIBS = -lm

SRC_DIR = src
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=%.o)
TARGET = inflation_easy

# Optional OpenMP flags
ifeq ($(shell $(CXX) -fopenmp -dM -E - < /dev/null > /dev/null 2>&1 && echo OK),OK)
	CXXFLAGS += -fopenmp
	LIBS += -fopenmp
endif

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)
	@rm -f $(OBJS)

%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/parameters.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@echo "Cleaning build files..."
	@rm -f $(OBJS) $(TARGET)
