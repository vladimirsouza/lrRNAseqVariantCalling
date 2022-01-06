#!/bin/bash



### inputs
REF=/home/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta
BAM_DIR=/home/vbarbo/project_2021/datasets/wtc11/manipulate_data
DEEPVARIANT_DIR=/home/vbarbo/project_2021/datasets/wtc11/methods_to_comp/deepvariant
THREADS=30



############################################
###### call variants with DeepVariant ######
############################################


###
### DV (alone)
###

singularity exec --bind $DEEPVARIANT_DIR,/usr/lib/locale/ \
  /home/vbarbo/programs/deepvariant_singularity/deepvariant-1.1.0.simg \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type PACBIO \
  --ref $REF \
  --reads $BAM_DIR/aln_s.bam \
  --output_vcf $DEEPVARIANT_DIR/no_sncr/deepvariant_calls.vcf \
  --num_shards $THREADS

# filter only PASS
bcftools view -f PASS \
  $DEEPVARIANT_DIR/no_sncr/deepvariant_calls.vcf \
  > $DEEPVARIANT_DIR/no_sncr/deepvariant_calls_pass.vcf

# compressing and indexing
bgzip -c $DEEPVARIANT_DIR/no_sncr/deepvariant_calls_pass.vcf \
  > $DEEPVARIANT_DIR/no_sncr/deepvariant_calls_pass.vcf.gz
bcftools index $DEEPVARIANT_DIR/no_sncr/deepvariant_calls_pass.vcf.gz
tabix -p vcf $DEEPVARIANT_DIR/no_sncr/deepvariant_calls_pass.vcf.gz




###
### SNCR + DV
###

singularity exec --bind $DEEPVARIANT_DIR,/usr/lib/locale/ \
  /home/vbarbo/programs/deepvariant_singularity/deepvariant-1.1.0.simg \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type PACBIO \
  --ref $REF \
  --reads $BAM_DIR/aln_sncr.bam \
  --output_vcf $DEEPVARIANT_DIR/yes_sncr_no_fc/deepvariant_calls.vcf \
  --num_shards $THREADS

# filter only PASS
bcftools view -f PASS \
  $DEEPVARIANT_DIR/yes_sncr_no_fc/deepvariant_calls.vcf \
  > $DEEPVARIANT_DIR/yes_sncr_no_fc/deepvariant_calls_pass.vcf

# compressing and indexing
bgzip -c $DEEPVARIANT_DIR/yes_sncr_no_fc/deepvariant_calls_pass.vcf \
  > $DEEPVARIANT_DIR/yes_sncr_no_fc/deepvariant_calls_pass.vcf.gz
bcftools index $DEEPVARIANT_DIR/yes_sncr_no_fc/deepvariant_calls_pass.vcf.gz
tabix -p vcf $DEEPVARIANT_DIR/yes_sncr_no_fc/deepvariant_calls_pass.vcf.gz





###
### SNCR + FC + DV
###

singularity exec --bind $DEEPVARIANT_DIR,/usr/lib/locale/ \
  /home/vbarbo/programs/deepvariant_singularity/deepvariant-1.1.0.simg \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type PACBIO \
  --ref $REF \
  --reads $BAM_DIR/aln_sncr_fc.bam \
  --output_vcf $DEEPVARIANT_DIR/deepvariant_calls.vcf \
  --num_shards $THREADS

# filter only PASS
bcftools view -f PASS \
  $DEEPVARIANT_DIR/deepvariant_calls.vcf \
  > $DEEPVARIANT_DIR/deepvariant_calls_pass.vcf

# compressing and indexing
bgzip -c $DEEPVARIANT_DIR/deepvariant_calls_pass.vcf \
  > $DEEPVARIANT_DIR/deepvariant_calls_pass.vcf.gz
bcftools index $DEEPVARIANT_DIR/deepvariant_calls_pass.vcf.gz
tabix -p vcf $DEEPVARIANT_DIR/deepvariant_calls_pass.vcf.gz


