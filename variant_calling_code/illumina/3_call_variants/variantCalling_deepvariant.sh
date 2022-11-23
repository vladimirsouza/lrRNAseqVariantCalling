#!/bin/bash


####### don't add read groups and don't mark/remove duplicates



### inputs
REF=/home/vbarbo/project_2021/paper_analysis/reference/genome/GRCh38.p13_all_chr.fasta
OUTPUT_DIR=/home/vbarbo/project_2021/paper_analysis/secondary_analyses/short_reads/files/variant_calling_from_short_reads/deepvariant
THREADS=40



############################################
###### call variants with DeepVariant ######
############################################


# DeepVariant alone (without BAM manipulation with SplitNCigarReads+flagCorrection) was killed after about 3 weeks running.


###
### DV (alone)
###
INPUT_BAM=/home/vbarbo/project_2021/paper_analysis/wtc11/ase_analysis/star_aln/GSE175048_Aligned.sortedByCoord.out_primary_readGroupAddedSameOrder.bam
mkdir -p ${OUTPUT_DIR}/dv/intermediate_results_dir
singularity exec --bind /usr/lib/locale/:/usr/lib/locale/ \
  /home/vbarbo/programs/deepvariant_singularity/deepvariant-1.1.0.simg \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type WES \
  --ref ${REF} \
  --reads ${INPUT_BAM} \
  --output_vcf ${OUTPUT_DIR}/dv/deepvariant_calls.vcf \
  --num_shards ${THREADS} \
  --intermediate_results_dir ${OUTPUT_DIR}/dv/intermediate_results_dir
rm -r ${OUTPUT_DIR}/dv/intermediate_results_dir

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
INPUT_BAM=/home/vbarbo/project_2021/paper_analysis/secondary_analyses/short_reads/files/data_manipulation/GSE175048_Aligned.sortedByCoord.out_primary_readGroupAddedSameOrder_sncr.bam
mkdir -p ${OUTPUT_DIR}/dv_sncr/intermediate_results_dir
singularity exec --bind /usr/lib/locale/:/usr/lib/locale/ \
  /home/vbarbo/programs/deepvariant_singularity/deepvariant-1.1.0.simg \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type WES \
  --ref ${REF} \
  --reads ${INPUT_BAM} \
  --output_vcf ${OUTPUT_DIR}/dv_sncr/deepvariant_calls.vcf \
  --num_shards ${THREADS} \
  --intermediate_results_dir ${OUTPUT_DIR}/dv_sncr/intermediate_results_dir
rm -r ${OUTPUT_DIR}/dv_sncr/intermediate_results_dir

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
INPUT_BAM=/home/vbarbo/project_2021/paper_analysis/secondary_analyses/short_reads/files/data_manipulation/GSE175048_Aligned.sortedByCoord.out_primary_readGroupAddedSameOrder_sncr_fc.bam
mkdir -p ${OUTPUT_DIR}/dv_sncr_fc/intermediate_results_dir
singularity exec --bind /usr/lib/locale/:/usr/lib/locale/ \
  /home/vbarbo/programs/deepvariant_singularity/deepvariant-1.1.0.simg \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type WES \
  --ref ${REF} \
  --reads ${INPUT_BAM} \
  --output_vcf ${OUTPUT_DIR}/dv_sncr_fc/deepvariant_calls.vcf \
  --num_shards ${THREADS} \
  --intermediate_results_dir ${OUTPUT_DIR}/dv_sncr_fc/intermediate_results_dir
rm -r ${OUTPUT_DIR}/dv_sncr_fc/intermediate_results_dir

# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls.vcf \
  > $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls_pass.vcf

# compressing and indexing
bgzip -c $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls_pass.vcf \
  > $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls_pass.vcf.gz
bcftools index $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/dv_sncr_fc/deepvariant_calls_pass.vcf.gz


