#!/bin/bash


####### don't add read groups and don't mark/remove duplicates



### inputs
REF=/home/vbarbo/project_2021/paper_analysis/reference/genome/GRCh38.p13_all_chr.fasta
BAM_DIR=/home/vbarbo/project_2021/paper_analysis/jurkat/data_manipulation
OUTPUT_DIR=/home/vbarbo/project_2021/paper_analysis/jurkat/variant_calling_from_isoseq/deepvariant
THREADS=30



############################################
###### call variants with DeepVariant ######
############################################


###
### DV (alone)
###

singularity exec --bind $OUTPUT_DIR/dv/,/usr/lib/locale/ \
  /home/vbarbo/programs/deepvariant_singularity/deepvariant-1.1.0.simg \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type PACBIO \
  --ref $REF \
  --reads $BAM_DIR/aln_s.bam \
  --output_vcf $OUTPUT_DIR/dv/deepvariant_calls.vcf \
  --num_shards $THREADS

# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/dv/deepvariant_calls.vcf \
  > $OUTPUT_DIR/dv/deepvariant_calls_pass.vcf

# compressing and indexing
bgzip -c $OUTPUT_DIR/dv/deepvariant_calls_pass.vcf \
  > $OUTPUT_DIR/dv/deepvariant_calls_pass.vcf.gz
bcftools index $OUTPUT_DIR/dv/deepvariant_calls_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/dv/deepvariant_calls_pass.vcf.gz




###
### SNCR + DV
###

singularity exec --bind $OUTPUT_DIR/dv_sncr/,/usr/lib/locale/ \
  /home/vbarbo/programs/deepvariant_singularity/deepvariant-1.1.0.simg \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type PACBIO \
  --ref $REF \
  --reads $BAM_DIR/aln_sncr.bam \
  --output_vcf $OUTPUT_DIR/dv_sncr/deepvariant_calls.vcf \
  --num_shards $THREADS

# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/dv_sncr/deepvariant_calls.vcf \
  > $OUTPUT_DIR/dv_sncr/deepvariant_calls_pass.vcf

# compressing and indexing
bgzip -c $OUTPUT_DIR/dv_sncr/deepvariant_calls_pass.vcf \
  > $OUTPUT_DIR/dv_sncr/deepvariant_calls_pass.vcf.gz
bcftools index $OUTPUT_DIR/dv_sncr/deepvariant_calls_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/dv_sncr/deepvariant_calls_pass.vcf.gz





###
### SNCR + FC + DV
###

singularity exec --bind $OUTPUT_DIR/dv_sncr_fc/,/usr/lib/locale/ \
  /home/vbarbo/programs/deepvariant_singularity/deepvariant-1.1.0.simg \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type PACBIO \
  --ref $REF \
  --reads $BAM_DIR/aln_sncr_fc.bam \
  --output_vcf $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls.vcf \
  --num_shards $THREADS

# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls.vcf \
  > $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls_pass.vcf

# compressing and indexing
bgzip -c $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls_pass.vcf \
  > $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls_pass.vcf.gz
bcftools index $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls_pass.vcf.gz


