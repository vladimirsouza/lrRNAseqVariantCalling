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
PATH_TO_REPO=/home/vbarbo/project_2021/projects/lrRNAseqVariantCalling




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

# in case the concatenate vcf ends up with two different variants at a same site,
# keeping the one that shows the highest QUAL value.
Rscript $PATH_TO_REPO/removeRepeatedLowerQualSites.r \
  $OUTPUT_DIR/mix/pileup_pass_mix.recode.vcf.gz \
  $OUTPUT_DIR/mix/pileup_pass_mix_norep.recode.vcf.gz



### delete temporary files
rm -r \
  $OUTPUT_DIR/c3/log/ \
  $OUTPUT_DIR/c3/tmp/ \
  $OUTPUT_DIR/c3_sncr/log/ \
  $OUTPUT_DIR/c3_sncr/tmp/ \
  $OUTPUT_DIR/c3_sncr_fc/log/ \
  $OUTPUT_DIR/c3_sncr_fc/tmp/

