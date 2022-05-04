CXX = clang++-brew

DBGFLAGS = \
  -g \
  -DNDEBUG

OPTFLAGS = \
  -O2 \
  -march=x86-64 \
  -mtune=znver1 \
  -mno-avx \
  -msse4

INCLUDES = \
  -isystem include

CXXFLAGS = \
  -std=c++17 \
  -pedantic \
  -Wall \
  -Wextra \
  -Wconversion \
  -Wsign-conversion \
  -Wshadow \
  -Wno-c99-extensions \
  $(DBGFLAGS) \
  $(OPTFLAGS) \
  $(INCLUDES)

LIBS = \
  -lhdf5 \
  -lz \
  -lprofiler

PRODUCTS = \
  simulation_spindle \
  simulation_interphase \
  simulation_fine_sampling

# Sources
COMMON_SOURCES = $(shell find src/simulation_common -name "*.cc")
COMMON_OBJECTS = $(COMMON_SOURCES:.cc=.o)

SPINDLE_SOURCES = $(shell find src/simulation_spindle -name "*.cc")
SPINDLE_OBJECTS = $(SPINDLE_SOURCES:.cc=.o)

INTERPHASE_SOURCES = $(shell find src/simulation_interphase -name "*.cc")
INTERPHASE_OBJECTS = $(INTERPHASE_SOURCES:.cc=.o)

FINE_SAMPLING_SOURCES = $(shell find src/simulation_fine_sampling -name "*.cc")
FINE_SAMPLING_OBJECTS = $(FINE_SAMPLING_SOURCES:.cc=.o)

SOURCES = $(COMMON_SOURCES) $(SPINDLE_SOURCES) $(INTERPHASE_SOURCES) $(FINE_SAMPLING_SOURCES)
OBJECTS = $(COMMON_OBJECTS) $(SPINDLE_OBJECTS) $(INTERPHASE_OBJECTS) $(FINE_SAMPLING_OBJECTS)
ARTIFACTS = $(PRODUCTS) $(OBJECTS)


.PHONY: all clean depends
.SUFFIXES: .cc

all: $(PRODUCTS)
	@:

clean:
	rm -f $(ARTIFACTS)

depends:
	for src in $(SOURCES); do \
	    $(CXX) $(CXXFLAGS) -MM -MF- -MT $${src%.cc}.o $$src; \
	done > depends.mk

.cc.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<

simulation_spindle: $(SPINDLE_OBJECTS) $(COMMON_OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

simulation_interphase: $(INTERPHASE_OBJECTS) $(COMMON_OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

simulation_fine_sampling: $(FINE_SAMPLING_OBJECTS) $(COMMON_OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

-include depends.mk
