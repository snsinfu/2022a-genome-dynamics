CPPFLAGS = -xc -std=c99 -P
PRODUCTS = defaults.py


.PHONY: all clean

all: $(PRODUCTS)
	@:

clean:
	rm -f $(PRODUCTS)


defaults.py: defaults.py.in
	$(CPP) $(CPPFLAGS) -o $@ $<

defaults.py: ../config_entries.inc
