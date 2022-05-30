CXXFLAGS = \
  $(WARNFLAGS) \
  $(OPTFLAGS) \
  $(DBGFLAGS) \
  $(INCLUDES)

WARNFLAGS = \
  -std=c++20 \
  -Wpedantic \
  -Wall \
  -Wextra \
  -Wconversion \
  -Wpedantic \
  -Wconversion \
  -Wsign-conversion \
  -Wshadow

OPTFLAGS = \
  -O2

DBGFLAGS = \
  -g

INCLUDES = \
  -isystem ../../submodules/github.com/danielaparker/jsoncons/include \
  -isystem ../../submodules/github.com/snsinfu/cxx-getopt \
  -isystem ../../submodules/github.com/snsinfu/h5/include \
  -isystem ../../submodules/github.com/snsinfu/micromd/include

LIBS = \
  -lhdf5 \
  -lprofiler

SOURCES = \
  $(wildcard *.cpp) \
  $(wildcard forces/*.cpp) \
  $(wildcard glues/*.cpp) \
  $(wildcard loops/*.cpp) \
  $(wildcard inits/*.cpp)

OBJECTS = \
  $(SOURCES:.cpp=.o)

PRODUCTS = \
  main

ARTIFACTS = \
  $(PRODUCTS) \
  $(OBJECTS) \
  _depends.mk


.PHONY: all clean depends

all: $(PRODUCTS)
	@:

clean:
	rm -f $(ARTIFACTS)

depends:
	for src in $(SOURCES); do \
	    $(CXX) $(CXXFLAGS) -c -MM -MF- -MT $${src%.*}.o $${src}; \
	done > _depends.mk

main: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

-include _depends.mk
-include _overrides.mk
