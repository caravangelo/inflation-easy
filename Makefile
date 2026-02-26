CXX      = c++
CXXFLAGS = -std=c++17 -O3 -Wall
LDFLAGS  =
LIBS     = -lm

SRC_DIR = src
SRCS    = $(wildcard $(SRC_DIR)/*.cpp)
OBJS    = $(SRCS:$(SRC_DIR)/%.cpp=%.o)
TARGET  = inflation_easy

# ---------- macOS exception: prefer Homebrew LLVM for OpenMP ----------
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  LLVM_PREFIX := $(shell brew --prefix llvm 2>/dev/null)
  ifneq ($(LLVM_PREFIX),)
    # Use brew's clang++ and wire OpenMP includes/libs/rpath
    CXX      := $(LLVM_PREFIX)/bin/clang++
    CXXFLAGS += -fopenmp -I$(LLVM_PREFIX)/include
    LDFLAGS  += -L$(LLVM_PREFIX)/lib -Wl,-rpath,$(LLVM_PREFIX)/lib
    # libomp is pulled in by -fopenmp with clang; no extra LIBS needed
  else
    # Apple clang doesn't support -fopenmp; leave OpenMP off
  endif
else
  # ---------- Non-macOS: your original OpenMP probe ----------
  ifeq ($(shell $(CXX) -fopenmp -dM -E - < /dev/null > /dev/null 2>&1 && echo OK),OK)
    CXXFLAGS += -fopenmp
    LIBS     += -fopenmp
  endif
endif

.PHONY: all clean dev-regression-main-n16

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)
	@rm -f $(OBJS)

# Minimal change: add ffteasy.hpp so objects rebuild when the header changes
%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/parameters.h $(SRC_DIR)/ffteasy.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@echo "Cleaning build files..."
	@rm -f $(OBJS) $(TARGET)

# Maintainer-only regression check: compare HEAD vs main at N=16.
dev-regression-main-n16:
	python3 tests/regression_main_n16.py --repo . --main-ref main --params params.numerical.txt
