#!/bin/bash


####### don't add read groups and don't mark/remove duplicates


### inputs
REF=/home/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta
OUTPUT_ROOT_DIR=/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls
OUTPUT_DIR=/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/noMarkDuplicate
THREADS=30



############################################
###### call variants with DeepVariant ######
############################################


###
### DV (alone)
###

singularity exec --bind $OUTPUT_DIR,/usr/lib/locale/ \
  /home/vbarbo/programs/deepvariant_singularity/deepvariant-1.1.0.simg \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type PACBIO \
  --ref $REF \
  --reads $OUTPUT_ROOT_DIR/aln.bam \
  --output_vcf $OUTPUT_DIR/dv_alone/deepvariant_calls.vcf \
  --num_shards $THREADS

# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/dv_alone/deepvariant_calls.vcf \
  > $OUTPUT_DIR/dv_alone/deepvariant_calls_pass.vcf

# compressing and indexing
bgzip -c $OUTPUT_DIR/dv_alone/deepvariant_calls_pass.vcf \
  > $OUTPUT_DIR/dv_alone/deepvariant_calls_pass.vcf.gz
bcftools index $OUTPUT_DIR/dv_alone/deepvariant_calls_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/dv_alone/deepvariant_calls_pass.vcf.gz




###
### SNCR + DV
###

singularity exec --bind $OUTPUT_DIR,/usr/lib/locale/ \
  /home/vbarbo/programs/deepvariant_singularity/deepvariant-1.1.0.simg \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type PACBIO \
  --ref $REF \
  --reads $OUTPUT_DIR/aln_split.bam \
  --output_vcf $OUTPUT_DIR/no_flagCorrection/deepvariant_calls.vcf \
  --num_shards $THREADS

# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/no_flagCorrection/deepvariant_calls.vcf \
  > $OUTPUT_DIR/no_flagCorrection/deepvariant_calls_pass.vcf

# compressing and indexing
bgzip -c $OUTPUT_DIR/no_flagCorrection/deepvariant_calls_pass.vcf \
  > $OUTPUT_DIR/no_flagCorrection/deepvariant_calls_pass.vcf.gz
bcftools index $OUTPUT_DIR/no_flagCorrection/deepvariant_calls_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/no_flagCorrection/deepvariant_calls_pass.vcf.gz





###
### SNCR + FC + DV
###

singularity exec --bind $OUTPUT_DIR,/usr/lib/locale/ \
  /home/vbarbo/programs/deepvariant_singularity/deepvariant-1.1.0.simg \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type PACBIO \
  --ref $REF \
  --reads $OUTPUT_DIR/aln_split_flagCorrection.bam \
  --output_vcf $OUTPUT_DIR/deepvariant_calls.vcf \
  --num_shards $THREADS

# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/deepvariant_calls.vcf \
  > $OUTPUT_DIR/deepvariant_calls_pass.vcf

# compressing and indexing
bgzip -c $OUTPUT_DIR/deepvariant_calls_pass.vcf \
  > $OUTPUT_DIR/deepvariant_calls_pass.vcf.gz
bcftools index $OUTPUT_DIR/deepvariant_calls_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/deepvariant_calls_pass.vcf.gz

