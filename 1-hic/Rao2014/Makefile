# Contact maps from 2014-Rao in-situ Hi-C experiments.
#
# NOTE: We don't use HMEC, HUVEC, K562, KBM7 and NHEK as these cell lines have
# irregular karyotypes, which makes analysis difficult.

DATA_DIRS = \
  data/raw \
  data/processed

# Raw *.hic files from GEO
RAW_FILES = \
  GSE63525_CH12-LX_combined_30.hic \
  GSE63525_GM12878_dilution_combined_30.hic \
  GSE63525_GM12878_insitu_primary+replicate_combined_30.hic \
  GSE63525_GM12878_insitu_primary+replicate_combined_30.hic \
  GSE63525_GM12878_diploid_maternal.hic \
  GSE63525_GM12878_diploid_paternal.hic \
  GSE63525_GM12878_insitu_replicate_30.hic \
  GSE63525_IMR90_combined_30.hic

# Example: GSE63525_IMR90_combined_30.mcool
COOL_FILES = \
  $(addsuffix .mcool, $(basename $(RAW_FILES)))

TARGETS = \
  $(addprefix data/interim/, $(COOL_FILES)) \
  data/raw/mm9.chrom.sizes \
  data/raw/hg19.chrom.sizes


.PHONY: all clean

all: $(TARGETS)
	@:

clean:
	rm -f $(TARGETS)

data/raw/GS%:
	curl -fsSLo "$@" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=$(shell echo '$@' | egrep -o 'GS[EM][0-9]+')&format=file&file=$(shell echo '$@' | sed 's/+/%2b/g')"

data/raw/hg19.chrom.sizes:
	curl -fsSLo "$@" http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes

data/raw/mm9.chrom.sizes:
	curl -fsSLo "$@" http://hgdownload.cse.ucsc.edu/goldenpath/mm9/bigZips/mm9.chrom.sizes

data/interim/%.mcool: data/raw/%.hic
	hic2cool convert "$^" "$@"
