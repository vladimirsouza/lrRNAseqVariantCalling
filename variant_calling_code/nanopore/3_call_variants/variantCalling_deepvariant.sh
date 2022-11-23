#!/bin/bash


####### don't add read groups and don't mark/remove duplicates



### inputs
REF=/home/vbarbo/project_2021/paper_analysis/reference/genome/GRCh38.p13_all_chr.fasta
BAM_DIR=/home/vbarbo/project_2021/paper_analysis/secondary_analyses/nanopore/files/ENCFF961HLO/data_manipulation
OUTPUT_DIR=/home/vbarbo/project_2021/paper_analysis/secondary_analyses/nanopore/files/ENCFF961HLO/variant_calling_from_nanopore/deepvariant
THREADS=10



############################################
###### call variants with DeepVariant ######
############################################


###
### DV (alone)
###

mkdir $OUTPUT_DIR/dv/intermediate_results_dir
singularity exec --bind /usr/lib/locale/:/usr/lib/locale/ \
  /home/vbarbo/programs/deepvariant_singularity/deepvariant-1.1.0.simg \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type WES \
  --ref $REF \
  --reads $BAM_DIR/aln_s.bam \
  --output_vcf $OUTPUT_DIR/dv/deepvariant_calls.vcf \
  --num_shards $THREADS \
  --intermediate_results_dir $OUTPUT_DIR/dv/intermediate_results_dir
rm -r $OUTPUT_DIR/dv/intermediate_results_dir


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
mkdir $OUTPUT_DIR/dv_sncr/intermediate_results_dir
singularity exec --bind $OUTPUT_DIR/dv_sncr/,/usr/lib/locale/ \
  /home/vbarbo/programs/deepvariant_singularity/deepvariant-1.1.0.simg \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type WES \
  --ref $REF \
  --reads $BAM_DIR/aln_sncr.bam \
  --output_vcf $OUTPUT_DIR/dv_sncr/deepvariant_calls.vcf \
  --num_shards $THREADS \
  --intermediate_results_dir $OUTPUT_DIR/dv_sncr/intermediate_results_dir
rm -r $OUTPUT_DIR/dv_sncr/intermediate_results_dir

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
mkdir $OUTPUT_DIR/dv_sncr_fc/intermediate_results_dir
singularity exec --bind $OUTPUT_DIR/dv_sncr_fc/,/usr/lib/locale/ \
  /home/vbarbo/programs/deepvariant_singularity/deepvariant-1.1.0.simg \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type WES \
  --ref $REF \
  --reads $BAM_DIR/aln_sncr_fc.bam \
  --output_vcf $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls.vcf \
  --num_shards $THREADS \
  --intermediate_results_dir $OUTPUT_DIR/dv_sncr_fc/intermediate_results_dir
rm -r $OUTPUT_DIR/dv_sncr_fc/intermediate_results_dir

# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls.vcf \
  > $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls_pass.vcf

# compressing and indexing
bgzip -c $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls_pass.vcf \
  > $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls_pass.vcf.gz
bcftools index $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls_pass.vcf.gz


