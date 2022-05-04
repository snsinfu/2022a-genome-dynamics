GENERATED = \
  Chains-PureAB-L20_N50.tsv \
  Chains-PureAB-L20_N200.tsv \
  Chains-SemiAB-L20_N200.tsv \
  Chains-PureAB-L10_N400.tsv \
  Chains-SemiAB-L10_N400.tsv \
  Chains-PureAB-L50_N80.tsv \
  Particles-N4000.tsv

TARGETS = \
  $(addprefix data/generated/, $(GENERATED))


.PHONY: all clean

all: $(TARGETS)
	@:

clean:
	rm -f $(TARGETS_GENERATED)

data/generated/Chains-PureAB-L20_N50.tsv:
	scripts/make_chain_definition pure 20 50 > $@

data/generated/Chains-PureAB-L20_N200.tsv:
	scripts/make_chain_definition pure 20 200 > $@

data/generated/Chains-SemiAB-L20_N200.tsv:
	scripts/make_chain_definition semi 20 200 > $@

data/generated/Chains-PureAB-L10_N400.tsv:
	scripts/make_chain_definition pure 10 400 > $@

data/generated/Chains-SemiAB-L10_N400.tsv:
	scripts/make_chain_definition semi 10 400 > $@

data/generated/Chains-PureAB-L50_N80.tsv:
	scripts/make_chain_definition pure 50 80 > $@

data/generated/Particles-N4000.tsv:
	scripts/make_particle_definition 4000 > $@
