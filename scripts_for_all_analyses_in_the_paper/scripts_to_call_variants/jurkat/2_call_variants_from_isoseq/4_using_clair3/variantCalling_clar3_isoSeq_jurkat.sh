#!/bin/bash



###################################
###### run clair3 on iso-seq ######
###################################

conda deactivate
conda activate clair3


REF="/home/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta"
THREADS=10


### Clair3 (alone)
INPUT_BAM="/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/aln.bam"
OUTPUT_DIR="/home/vbarbo/project_2021/datasets/gloria_data/analysis/clair3_isoSeq/alone"
# call variants
/home/vbarbo/programs/Clair3/run_clair3.sh \
  --bam_fn=$INPUT_BAM \
  --ref_fn=$REF \
  --threads=$THREADS \
  --platform="hifi" \
  --model_path="/home/vbarbo/programs/Clair3/models/hifi" \
  --output=$OUTPUT_DIR
# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/pileup.vcf.gz \
  > $OUTPUT_DIR/pileup_pass.vcf
# compress and index
bgzip -c $OUTPUT_DIR/pileup_pass.vcf \
  > $OUTPUT_DIR/pileup_pass.vcf.gz
bcftools index $OUTPUT_DIR/pileup_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/pileup_pass.vcf.gz
rm $OUTPUT_DIR/pileup_pass.vcf


### SNCR + Clair3
INPUT_BAM="/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/noMarkDuplicate/aln_split.bam"
OUTPUT_DIR="/home/vbarbo/project_2021/datasets/gloria_data/analysis/clair3_isoSeq/sncr"
/home/vbarbo/programs/Clair3/run_clair3.sh \
  --bam_fn=$INPUT_BAM \
  --ref_fn=$REF \
  --threads=$THREADS \
  --platform="hifi" \
  --model_path="/home/vbarbo/programs/Clair3/models/hifi" \
  --output=$OUTPUT_DIR
# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/pileup.vcf.gz \
  > $OUTPUT_DIR/pileup_pass.vcf
# compress and index
bgzip -c $OUTPUT_DIR/pileup_pass.vcf \
  > $OUTPUT_DIR/pileup_pass.vcf.gz
bcftools index $OUTPUT_DIR/pileup_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/pileup_pass.vcf.gz
rm $OUTPUT_DIR/pileup_pass.vcf


### SNCR + FC + Clair3
INPUT_BAM="/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/noMarkDuplicate/aln_split_flagCorrection.bam"
OUTPUT_DIR="/home/vbarbo/project_2021/datasets/gloria_data/analysis/clair3_isoSeq/sncr_fc"
/home/vbarbo/programs/Clair3/run_clair3.sh \
  --bam_fn=$INPUT_BAM \
  --ref_fn=$REF \
  --threads=$THREADS \
  --platform="hifi" \
  --model_path="/home/vbarbo/programs/Clair3/models/hifi" \
  --output=$OUTPUT_DIR
# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/pileup.vcf.gz \
  > $OUTPUT_DIR/pileup_pass.vcf
# compress and index
bgzip -c $OUTPUT_DIR/pileup_pass.vcf \
  > $OUTPUT_DIR/pileup_pass.vcf.gz
bcftools index $OUTPUT_DIR/pileup_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/pileup_pass.vcf.gz
rm $OUTPUT_DIR/pileup_pass.vcf






### SNPs from Clair3 (alone), indels from SNCR + FC+ Clair3
OUTPUT_DIR="/home/vbarbo/project_2021/datasets/gloria_data/analysis/clair3_isoSeq/mix"
VCF_ALONE=/home/vbarbo/project_2021/datasets/gloria_data/analysis/clair3_isoSeq/alone/pileup_pass.vcf.gz
VCF_SNCR_FC=/home/vbarbo/project_2021/datasets/gloria_data/analysis/clair3_isoSeq/sncr_fc/pileup_pass.vcf.gz


mkdir $OUTPUT_DIR

# snps
vcftools --gzvcf $VCF_ALONE \
  --out $OUTPUT_DIR/pileup_pass_snp \
  --remove-indels --recode --recode-INFO-all
bgzip $OUTPUT_DIR/pileup_pass_snp.recode.vcf
tabix -p vcf $OUTPUT_DIR/pileup_pass_snp.recode.vcf.gz

# indels
vcftools --gzvcf $VCF_SNCR_FC \
  --out $OUTPUT_DIR/pileup_pass_indel \
  --keep-only-indels --recode --recode-INFO-all
bgzip $OUTPUT_DIR/pileup_pass_indel.recode.vcf
tabix -p vcf $OUTPUT_DIR/pileup_pass_indel.recode.vcf.gz

# concatenate vcf files
bcftools concat \
  $OUTPUT_DIR/pileup_pass_snp.recode.vcf.gz \
  $OUTPUT_DIR/pileup_pass_indel.recode.vcf.gz \
  -o $OUTPUT_DIR/pileup_pass_mix.recode.vcf.gz \
  -O z -D -a

### to remove one of the (different) variants with same positions
### >>> This is done in R
library(vcfR)
vcf <- read.vcfR("/home/vbarbo/project_2021/datasets/gloria_data/analysis/clair3_isoSeq/mix/pileup_pass_mix.recode.vcf.gz")
pos <- paste(vcf@fix[,1], vcf@fix[,2])
k <- duplicated(pos)
vcf_nodup <- vcf[!k]
write.vcf(vcf_nodup, file="/home/vbarbo/project_2021/datasets/gloria_data/analysis/clair3_isoSeq/mix/pileup_pass_mix_nodup.recode.vcf.gz")
### <<< close R
# for some reason, bcftools can't read vcf files produced by vcfR::write.vcf. 
# decompressing and compressing it again can solve the problem!
gunzip /home/vbarbo/project_2021/datasets/gloria_data/analysis/clair3_isoSeq/mix/pileup_pass_mix_nodup.recode.vcf.gz
bgzip /home/vbarbo/project_2021/datasets/gloria_data/analysis/clair3_isoSeq/mix/pileup_pass_mix_nodup.recode.vcf


