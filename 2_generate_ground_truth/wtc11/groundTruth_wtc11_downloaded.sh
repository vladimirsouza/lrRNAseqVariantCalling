#!/bin/bash


# We filter the high quality VCF file downloaded from Allen Institute to generate our ground truth.


ORIGINAL_WTC11_VCF=/home/vbarbo/project_2021/paper_analysis/wtc11/ground_truth/3546dc62_AH77TTBBXX_DS-229105_GCCAAT_recalibrated.vcf.gz
OUTPUT_DIR=/home/vbarbo/project_2021/paper_analysis/wtc11/ground_truth


cd $OUTPUT_DIR

### index
bcftools index \
  $ORIGINAL_WTC11_VCF

### subset only chromosomes chr1, chr2, ..., chr22, chrX, chrY
bcftools view \
  $ORIGINAL_WTC11_VCF \
  --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
  > $OUTPUT_DIR/3546dc62_AH77TTBBXX_DS-229105_GCCAAT_recalibrated_subsetChromosomes.vcf

### compress
bgzip -c $OUTPUT_DIR/3546dc62_AH77TTBBXX_DS-229105_GCCAAT_recalibrated_subsetChromosomes.vcf \
  > $OUTPUT_DIR/3546dc62_AH77TTBBXX_DS-229105_GCCAAT_recalibrated_subsetChromosomes.vcf.gz
rm $OUTPUT_DIR/3546dc62_AH77TTBBXX_DS-229105_GCCAAT_recalibrated_subsetChromosomes.vcf


### filter out non-PASS variants
bcftools view -f PASS \
  $OUTPUT_DIR/3546dc62_AH77TTBBXX_DS-229105_GCCAAT_recalibrated_subsetChromosomes.vcf.gz \
  > $OUTPUT_DIR/3546dc62_AH77TTBBXX_DS-229105_GCCAAT_recalibrated_subsetChromosomes_pass.vcf

### compress and index
bgzip -c $OUTPUT_DIR/3546dc62_AH77TTBBXX_DS-229105_GCCAAT_recalibrated_subsetChromosomes_pass.vcf \
  > $OUTPUT_DIR/3546dc62_AH77TTBBXX_DS-229105_GCCAAT_recalibrated_subsetChromosomes_pass.vcf.gz
rm $OUTPUT_DIR/3546dc62_AH77TTBBXX_DS-229105_GCCAAT_recalibrated_subsetChromosomes_pass.vcf
bcftools index $OUTPUT_DIR/3546dc62_AH77TTBBXX_DS-229105_GCCAAT_recalibrated_subsetChromosomes_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/3546dc62_AH77TTBBXX_DS-229105_GCCAAT_recalibrated_subsetChromosomes_pass.vcf.gz


