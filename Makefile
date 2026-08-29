CXX      = c++
CPPFLAGS ?=
CXXFLAGS = -std=c++17 -O3 -DNDEBUG -Wall
LDFLAGS  =
LIBS     = -lm
ENABLE_NATIVE ?= 0
ENABLE_LTO    ?= 0

SRC_DIR       = src
EVOLUTION_DIR = $(SRC_DIR)/evolution
SOURCE_DIRS   = $(SRC_DIR) $(EVOLUTION_DIR)
BUILD_DIR     = build
CONFIG_STAMP  = $(BUILD_DIR)/.build_config

SRCS   = $(foreach dir,$(SOURCE_DIRS),$(wildcard $(dir)/*.cpp))
OBJS   = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))
DEPS   = $(OBJS:.o=.d)
TARGET = inflation_easy

PROJECT_CPPFLAGS = -I$(SRC_DIR)

# ---------- macOS exception: prefer Homebrew LLVM for OpenMP ----------
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  LLVM_PREFIX := $(shell brew --prefix llvm 2>/dev/null)
  ifneq ($(LLVM_PREFIX),)
    LLVM_CXX := $(LLVM_PREFIX)/bin/clang++
    ifneq ($(wildcard $(LLVM_CXX)),)
      # Use Homebrew clang++ when it is installed, with OpenMP enabled.
      CXX      := $(LLVM_CXX)
      CXXFLAGS += -fopenmp -I$(LLVM_PREFIX)/include
      LDFLAGS  += -L$(LLVM_PREFIX)/lib -Wl,-rpath,$(LLVM_PREFIX)/lib
    endif
  endif
else
  # ---------- Non-macOS: enable OpenMP when the compiler supports it ----------
  ifeq ($(shell $(CXX) -fopenmp -dM -E - < /dev/null > /dev/null 2>&1 && echo OK),OK)
    CXXFLAGS += -fopenmp
    LIBS     += -fopenmp
  endif
endif

ifeq ($(ENABLE_NATIVE),1)
  CXXFLAGS += -march=native -mtune=native
endif

ifeq ($(ENABLE_LTO),1)
  CXXFLAGS += -flto
  LDFLAGS  += -flto
endif

.PHONY: all clean test test-spatial test-smoke test-sanitizers test-release dev-regression-main-n16 FORCE

all: $(TARGET)

# Rebuild when compiler or linker options change between make invocations.
# The timestamp changes only when the effective configuration changes.
FORCE:

$(CONFIG_STAMP): FORCE
	@mkdir -p $(dir $@)
	@printf '%s\n' \
	  'CXX=$(CXX)' \
	  'PROJECT_CPPFLAGS=$(PROJECT_CPPFLAGS)' \
	  'CPPFLAGS=$(CPPFLAGS)' \
	  'CXXFLAGS=$(CXXFLAGS)' \
	  'LDFLAGS=$(LDFLAGS)' \
	  'LIBS=$(LIBS)' > $@.tmp
	@if ! cmp -s $@.tmp $@; then mv $@.tmp $@; else rm -f $@.tmp; fi

$(OBJS): $(CONFIG_STAMP) Makefile

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

# Mirror the source tree under build/ and let the compiler track exact header
# dependencies for both the top-level sources and the evolution subsystem.
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(PROJECT_CPPFLAGS) $(CPPFLAGS) $(CXXFLAGS) \
	  -MMD -MP -MF $(@:.o=.d) -c $< -o $@

-include $(DEPS)

clean:
	@echo "Cleaning build files..."
	@rm -rf ./build
	@rm -f ./inflation_easy

# Verify the real-space operators and their Fourier eigenvalues for every
# supported spatial order. The temporary executables are kept outside the tree.
test-spatial:
	@for order in 2 4 6; do \
	  binary=$$(mktemp "/tmp/inflationeasy_spatial_test_$${order}.XXXXXX") || exit 1; \
	  if ! $(CXX) $(PROJECT_CPPFLAGS) $(CPPFLAGS) $(CXXFLAGS) \
	    -DSPATIAL_STENCIL_ORDER=$$order \
	    tests/spatial_discretization_test.cpp -o $$binary $(LDFLAGS) $(LIBS); then \
	    rm -f $$binary; exit 1; \
	  fi; \
	  if ! $$binary; then rm -f $$binary; exit 1; fi; \
	  rm -f $$binary; \
	done

# Fast checks suitable for every push and pull request.
test: test-spatial test-smoke

test-smoke:
	python3 tests/release_smoke.py --repo . --tier ci

# Focused memory/undefined-behaviour checks for historically delicate outputs.
test-sanitizers:
	python3 tests/release_smoke.py --repo . --tier sanitizers

# Broader pre-tag matrix. This is intentionally separate from the fast CI tier.
test-release:
	python3 tests/release_smoke.py --repo . --tier release

# Compatibility alias for the former branch-vs-main check. The old comparison
# becomes vacuous on main; use the strict clean-tree smoke suite instead.
dev-regression-main-n16:
	@echo "dev-regression-main-n16 is superseded by test-smoke"
	$(MAKE) test-smoke
