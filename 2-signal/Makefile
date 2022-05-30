MCOOLS = \
  GSE63525_GM12878_insitu_primary+replicate_combined_30.mcool \
  GSE63525_GM12878_diploid_maternal.mcool \
  GSE63525_GM12878_diploid_paternal.mcool \
  GSE63525_IMR90_combined_30.mcool \
  GSE63525_CH12-LX_combined_30.mcool \
  GSE104333_Rao-2017-untreated_combined_30.mcool \
  GSE104333_Rao-2017-treated_6hr_combined_30.mcool \
  GSE104333_Rao-2017-untreated_synchronized_combined_30.mcool \
  GSE104333_Rao-2017-treated_6hr_synchronized_combined_30.mcool

# Example: GSE63525_IMR90_combined_30
NAMES = \
  $(basename $(MCOOLS))

# Example: GSE63525_IMR90_combined_30-LIS-50kb.tsv
RAW_INTERACTION_TSVS = \
  $(addsuffix -LIS-50kb.tsv, $(NAMES)) \
  $(addsuffix -LIS-100kb.tsv, $(NAMES))

# Example: GSE63525_IMR90_combined_30-Alpha-10kb.tsv
RAW_ALPHA_TSVS = \
  $(addsuffix -Alpha-10kb.tsv, $(NAMES)) \
  $(addsuffix -Alpha-50kb.tsv, $(NAMES)) \
  $(addsuffix -Alpha-100kb.tsv, $(NAMES))

# Example: GSE63525_IMR90_combined_30-LIS-smooth_10x50kb-100kb.tsv
SMOOTH_INTERACTION_TSVS = \
  $(addsuffix -LIS-smooth_2x50kb-100kb.tsv, $(NAMES)) \
  $(addsuffix -LIS-smooth_10x50kb-100kb.tsv, $(NAMES)) \
  $(addsuffix -LIS-smooth_20x50kb-100kb.tsv, $(NAMES)) \
  $(addsuffix -LIS-smooth_50x50kb-100kb.tsv, $(NAMES))

# Genome models
GENOME_MODEL_TSVS = \
  Rao14-IMR90-Genome-rough.tsv \
  Rao14-IMR90-Genome-thres3_rough.tsv \
  Rao14-IMR90-Genome-smooth.tsv \
  Rao14-IMR90-Genome-thres0_smooth.tsv \
  Rao14-IMR90-Genome-thres2_smooth.tsv \
  Rao14-IMR90-Genome-thres3_smooth.tsv \
  Rao14-IMR90-Genome-thres4_smooth.tsv \
  Rao14-IMR90_hypernor-Genome-thres3_smooth.tsv \
  Rao14-IMR90_random-Genome.tsv \
  Rao14-IMR90_completely_random-Genome.tsv \
  Rao14-GM12878-Genome-rough.tsv \
  Rao14-GM12878-Genome-thres3_rough.tsv \
  Rao14-GM12878-Genome-smooth.tsv \
  Rao14-GM12878-Genome-thres0_smooth.tsv \
  Rao14-GM12878-Genome-thres2_smooth.tsv \
  Rao14-GM12878-Genome-thres3_smooth.tsv \
  Rao14-GM12878-Genome-thres4_smooth.tsv \
  Rao14-GM12878_hypernor-Genome-thres3_smooth.tsv \
  Rao14-GM12878_random-Genome.tsv \
  Rao17-HCT116-Genome-rough.tsv \
  Rao17-HCT116-Genome-smooth.tsv \
  Rao17-HCT116-Genome-thres3_smooth.tsv \
  Rao17-HCT116_treated-Genome-rough.tsv \
  Rao17-HCT116_treated-Genome-smooth.tsv \
  Rao17-HCT116_treated-Genome-thres3_smooth.tsv \
  Rao14-CH12LX-Genome.tsv

TARGETS = \
  $(addprefix data/interim/, $(RAW_INTERACTION_TSVS)) \
  $(addprefix data/interim/, $(RAW_ALPHA_TSVS)) \
  $(addprefix data/interim/, $(SMOOTH_INTERACTION_TSVS)) \
  $(addprefix data/processed/, $(GENOME_MODEL_TSVS))


.PHONY: all clean

all: $(TARGETS)
	@:

clean:
	rm -f $(TARGETS)

# LIS tracks

data/interim/%-LIS-25kb.tsv: data/external/%.mcool
	scripts/compute_interactions -b 25000 -o $@ $<

data/interim/%-LIS-50kb.tsv: data/external/%.mcool
	scripts/compute_interactions -b 50000 -o $@ $<

data/interim/%-LIS-100kb.tsv: data/external/%.mcool
	scripts/compute_interactions -b 100000 -o $@ $<

data/interim/%-LIS-smooth_2x50kb-100kb.tsv: data/interim/%-LIS-50kb.tsv
	scripts/downsample --rate 2 --window 2 -o $@ $<

data/interim/%-LIS-smooth_10x50kb-100kb.tsv: data/interim/%-LIS-50kb.tsv
	scripts/downsample --rate 2 --window 10 -o $@ $<

data/interim/%-LIS-smooth_20x50kb-100kb.tsv: data/interim/%-LIS-50kb.tsv
	scripts/downsample --rate 2 --window 20 -o $@ $<

data/interim/%-LIS-smooth_50x50kb-100kb.tsv: data/interim/%-LIS-50kb.tsv
	scripts/downsample --rate 2 --window 50 -o $@ $<

# Alpha tracks

data/interim/%-Alpha-10kb.tsv: data/external/%.mcool
	scripts/compute_local_alpha -w 50 -b 10000 -o $@ $<
	sed -i.bak "2,$$ s/^/chr/" $@
	rm $@.bak

data/interim/%-Alpha-50kb.tsv: data/external/%.mcool
	scripts/compute_local_alpha -w 50 -b 50000 -o $@ $<
	sed -i.bak "2,$$ s/^/chr/" $@
	rm $@.bak

data/interim/%-Alpha-100kb.tsv: data/external/%.mcool
	scripts/compute_local_alpha -w 50 -b 100000 -o $@ $<
	sed -i.bak "2,$$ s/^/chr/" $@
	rm $@.bak

# Genome models

data/processed/Rao14-IMR90-Genome-rough.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE63525_IMR90_combined_30-LIS-smooth_2x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --Xi -o $@ $^

data/processed/Rao14-IMR90-Genome-thres3_rough.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE63525_IMR90_combined_30-LIS-smooth_2x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --Xi --threshold 0.3 -o $@ $^

data/processed/Rao14-IMR90-Genome-smooth.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE63525_IMR90_combined_30-LIS-smooth_10x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --Xi -o $@ $^

data/processed/Rao14-IMR90-Genome-thres0_smooth.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE63525_IMR90_combined_30-LIS-smooth_10x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --Xi --threshold 0 -o $@ $^

data/processed/Rao14-IMR90-Genome-thres2_smooth.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE63525_IMR90_combined_30-LIS-smooth_10x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --Xi --threshold 0.2 -o $@ $^

data/processed/Rao14-IMR90-Genome-thres3_smooth.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE63525_IMR90_combined_30-LIS-smooth_10x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --Xi --threshold 0.3 -o $@ $^

data/processed/Rao14-IMR90-Genome-thres4_smooth.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE63525_IMR90_combined_30-LIS-smooth_10x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --Xi --threshold 0.4 -o $@ $^

data/processed/Rao14-IMR90_hypernor-Genome-thres3_smooth.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE63525_IMR90_combined_30-LIS-smooth_10x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --Xi --threshold 0.3 --hyperactive-nor -o $@ $^

data/processed/Rao14-IMR90_random-Genome.tsv: \
    data/processed/Rao14-IMR90-Genome-thres3_smooth.tsv
	scripts/randomize --seed 0 --preserve-structure -o $@ $^

data/processed/Rao14-IMR90_completely_random-Genome.tsv: \
    data/processed/Rao14-IMR90-Genome-thres3_smooth.tsv
	scripts/randomize --seed 0 --completely-random -o $@ $^

data/processed/Rao14-GM12878-Genome-rough.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE63525_GM12878_insitu_primary+replicate_combined_30-LIS-smooth_2x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --Xi -o $@ $^

data/processed/Rao14-GM12878-Genome-thres3_rough.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE63525_GM12878_insitu_primary+replicate_combined_30-LIS-smooth_2x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --Xi --threshold 0.3 -o $@ $^

data/processed/Rao14-GM12878-Genome-smooth.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE63525_GM12878_insitu_primary+replicate_combined_30-LIS-smooth_10x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --Xi -o $@ $^

data/processed/Rao14-GM12878-Genome-thres0_smooth.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE63525_GM12878_insitu_primary+replicate_combined_30-LIS-smooth_10x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --Xi --threshold 0 -o $@ $^

data/processed/Rao14-GM12878-Genome-thres2_smooth.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE63525_GM12878_insitu_primary+replicate_combined_30-LIS-smooth_10x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --Xi --threshold 0.2 -o $@ $^

data/processed/Rao14-GM12878-Genome-thres3_smooth.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE63525_GM12878_insitu_primary+replicate_combined_30-LIS-smooth_10x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --Xi --threshold 0.3 -o $@ $^

data/processed/Rao14-GM12878-Genome-thres4_smooth.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE63525_GM12878_insitu_primary+replicate_combined_30-LIS-smooth_10x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --Xi --threshold 0.4 -o $@ $^

data/processed/Rao14-GM12878_hypernor-Genome-thres3_smooth.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE63525_GM12878_insitu_primary+replicate_combined_30-LIS-smooth_10x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --Xi --threshold 0.3 --hyperactive-nor -o $@ $^

data/processed/Rao14-GM12878_random-Genome.tsv: \
    data/processed/Rao14-GM12878-Genome-thres3_smooth.tsv
	scripts/randomize --seed 0 --preserve-structure -o $@ $^

data/processed/Rao17-HCT116-Genome-rough.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE104333_Rao-2017-untreated_combined_30-LIS-smooth_2x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --criterion I1 -o $@ $^

data/processed/Rao17-HCT116-Genome-smooth.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE104333_Rao-2017-untreated_combined_30-LIS-smooth_10x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --criterion I1 -o $@ $^

data/processed/Rao17-HCT116-Genome-thres3_smooth.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE104333_Rao-2017-untreated_combined_30-LIS-smooth_10x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --threshold 0.3 -o $@ $^

data/processed/Rao17-HCT116_treated-Genome-rough.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE104333_Rao-2017-treated_6hr_combined_30-LIS-smooth_2x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --criterion I1 -o $@ $^

data/processed/Rao17-HCT116_treated-Genome-smooth.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE104333_Rao-2017-treated_6hr_combined_30-LIS-smooth_10x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --criterion I1 -o $@ $^

data/processed/Rao17-HCT116_treated-Genome-thres3_smooth.tsv: \
    data/external/genome/hg19-genome-100kb.tsv \
    data/interim/GSE104333_Rao-2017-treated_6hr_combined_30-LIS-smooth_10x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --threshold 0.3 -o $@ $^

data/processed/Rao14-CH12LX-Genome.tsv: \
    data/external/genome/mm9-genome-100kb.tsv \
    data/interim/GSE63525_CH12-LX_combined_30-LIS-smooth_10x50kb-100kb.tsv
	scripts/model_genome --diploid --Xa --Xi --threshold 0.3 -o $@ $^
