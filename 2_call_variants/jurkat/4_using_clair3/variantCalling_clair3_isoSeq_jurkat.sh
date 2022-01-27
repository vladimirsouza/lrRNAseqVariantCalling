#!/bin/bash



###################################
###### run clair3 on iso-seq ######
###################################


conda deactivate
conda activate clair3


REF=/home/vbarbo/project_2021/paper_analysis/reference/genome/GRCh38.p13_all_chr.fasta
INPUT_BAM_DIR=/home/vbarbo/project_2021/paper_analysis/jurkat/data_manipulation
OUTPUT_DIR=/home/vbarbo/project_2021/paper_analysis/jurkat/variant_calling_from_isoseq/clair3
THREADS=20




### Clair3 (alone)
# call variants
/home/vbarbo/programs/Clair3/run_clair3.sh \
  --bam_fn=$INPUT_BAM_DIR/aln.bam \
  --ref_fn=$REF \
  --threads=$THREADS \
  --platform="hifi" \
  --model_path="/home/vbarbo/programs/Clair3/models/hifi" \
  --output=$OUTPUT_DIR/c3
# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/c3/pileup.vcf.gz \
  > $OUTPUT_DIR/c3/pileup_pass.vcf
# compress and index
bgzip -c $OUTPUT_DIR/c3/pileup_pass.vcf \
  > $OUTPUT_DIR/c3/pileup_pass.vcf.gz
bcftools index $OUTPUT_DIR/c3/pileup_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/c3/pileup_pass.vcf.gz
rm $OUTPUT_DIR/c3/pileup_pass.vcf



### SNCR + Clair3
/home/vbarbo/programs/Clair3/run_clair3.sh \
  --bam_fn=$INPUT_BAM_DIR/aln_sncr.bam \
  --ref_fn=$REF \
  --threads=$THREADS \
  --platform="hifi" \
  --model_path="/home/vbarbo/programs/Clair3/models/hifi" \
  --output=$OUTPUT_DIR/c3_sncr
# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/c3_sncr/pileup.vcf.gz \
  > $OUTPUT_DIR/c3_sncr/pileup_pass.vcf
# compress and index
bgzip -c $OUTPUT_DIR/c3_sncr/pileup_pass.vcf \
  > $OUTPUT_DIR/c3_sncr/pileup_pass.vcf.gz
bcftools index $OUTPUT_DIR/c3_sncr/pileup_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/c3_sncr/pileup_pass.vcf.gz
rm $OUTPUT_DIR/c3_sncr/pileup_pass.vcf



### SNCR + FC + Clair3
/home/vbarbo/programs/Clair3/run_clair3.sh \
  --bam_fn=$INPUT_BAM_DIR/aln_sncr_fc.bam \
  --ref_fn=$REF \
  --threads=$THREADS \
  --platform="hifi" \
  --model_path="/home/vbarbo/programs/Clair3/models/hifi" \
  --output=$OUTPUT_DIR/c3_sncr_fc
# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/c3_sncr_fc/pileup.vcf.gz \
  > $OUTPUT_DIR/c3_sncr_fc/pileup_pass.vcf
# compress and index
bgzip -c $OUTPUT_DIR/c3_sncr_fc/pileup_pass.vcf \
  > $OUTPUT_DIR/c3_sncr_fc/pileup_pass.vcf.gz
bcftools index $OUTPUT_DIR/c3_sncr_fc/pileup_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/c3_sncr_fc/pileup_pass.vcf.gz
rm $OUTPUT_DIR/c3_sncr_fc/pileup_pass.vcf



### Clair3-mix:
### SNPs from Clair3 (alone), indels from SNCR + FC+ Clair3

# snps
vcftools --gzvcf $OUTPUT_DIR/c3/pileup_pass.vcf.gz \
  --out $OUTPUT_DIR/mix/pileup_pass_snp \
  --remove-indels --recode --recode-INFO-all
bgzip $OUTPUT_DIR/mix/pileup_pass_snp.recode.vcf
tabix -p vcf $OUTPUT_DIR/mix/pileup_pass_snp.recode.vcf.gz

# indels
vcftools --gzvcf $OUTPUT_DIR/c3_sncr_fc/pileup_pass.vcf.gz \
  --out $OUTPUT_DIR/mix/pileup_pass_indel \
  --keep-only-indels --recode --recode-INFO-all
bgzip $OUTPUT_DIR/mix/pileup_pass_indel.recode.vcf
tabix -p vcf $OUTPUT_DIR/mix/pileup_pass_indel.recode.vcf.gz

# concatenate vcf files
bcftools concat \
  $OUTPUT_DIR/mix/pileup_pass_snp.recode.vcf.gz \
  $OUTPUT_DIR/mix/pileup_pass_indel.recode.vcf.gz \
  -o $OUTPUT_DIR/mix/pileup_pass_mix.recode.vcf.gz \
  -O z -D -a

# to remove one of the (different) variants with same positions
Rscript -e "{
  library(vcfR)
  vcf <- read.vcfR('$OUTPUT_DIR/mix/pileup_pass_mix.recode.vcf.gz')
  pos <- paste(vcf@fix[,1], vcf@fix[,2])
  k <- duplicated(pos)
  vcf_nodup <- vcf[!k]
  write.vcf(vcf_nodup, file='$OUTPUT_DIR/mix/pileup_pass_mix_nodup.recode.vcf.gz')
}"

# for some reason, bcftools can't read vcf files produced by vcfR::write.vcf
# decompressing and compressing it again can solve the problem!!
gunzip $OUTPUT_DIR/mix/pileup_pass_mix_nodup.recode.vcf.gz
bgzip $OUTPUT_DIR/mix/pileup_pass_mix_nodup.recode.vcf


