INCLUDES = \
  -isystem ../../include

LIBS = \
  -lhdf5 \
  -lz

PRODUCT = rdf_analysis

SOURCES = $(wildcard *.cc)
OBJECTS = $(SOURCES:.cc=.o)

ARTIFACTS = \
  $(OBJECTS) \
  $(PRODUCT)


.PHONY: all clean

all: $(PRODUCT)
	@:

clean:
	rm -f $(ARTIFACTS)

$(PRODUCT): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

include ../Mk/common_cxx.mk
