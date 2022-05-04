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


.PHONY: depends
.SUFFIXES: .cc

.cc.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<

depends:
	for src in $(SOURCES); do \
	    $(CXX) $(CXXFLAGS) -MM -MF- -MT $${src%.cc}.o $$src; \
	done > depends.mk

-include depends.mk
