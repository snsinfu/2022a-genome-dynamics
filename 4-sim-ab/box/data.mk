GENERATED = \
  Chains-PureAB-L20_N50.tsv \
  Chains-PureAB-L20_N200.tsv \
  Chains-SemiAB-L20_N200.tsv \
  Chains-PureAB-L10_N400.tsv \
  Chains-SemiAB-L10_N400.tsv \
  Chains-PureAB-L20_Na150_Nb50.tsv \
  Chains-PureAB-La20_Na150_Lb10Nb150.tsv \
  Chains-PureAB-La20_Na100_Lb10Nb100.tsv

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

data/generated/Chains-PureAB-L20_Na150_Nb50.tsv:
	scripts/make_chain_definition_het 20 20 150 50 > $@

data/generated/Chains-PureAB-La20_Na150_Lb10Nb150.tsv:
	scripts/make_chain_definition_het 20 10 150 150 > $@

data/generated/Chains-PureAB-La20_Na100_Lb10Nb100.tsv:
	scripts/make_chain_definition_het 20 10 100 100 > $@
