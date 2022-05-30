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
  -isystem ../../include

override LDFLAGS += \
  -lz \
  -lhdf5

override CXXFLAGS += \
  -std=c++17 \
  -Wpedantic \
  -Wall \
  -Wextra \
  -Wconversion \
  -Wsign-conversion \
  -Wno-c99-extensions \
  $(DBGFLAGS) \
  $(OPTFLAGS) \
  $(INCLUDES)

PRODUCTS = \
  simulation

SOURCES = \
  $(wildcard *.cc)

OBJECTS = \
  $(SOURCES:.cc=.o)

ARTIFACTS = \
  $(PRODUCTS) \
  $(OBJECTS) \
  depends.mk


.PHONY: all clean depends

all: $(PRODUCTS)
	@:

clean:
	rm -f $(ARTIFACTS)

depends:
	for src in $(SOURCES); do \
	    $(CXX) $(CXXFLAGS) -MM -MF- -MT $${src%.cc}.o $$src; \
	done > depends.mk

simulation: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

-include depends.mk
