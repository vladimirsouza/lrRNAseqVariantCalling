#!/bin/bash



####################################
###### run clair3 on nanopore ######
####################################


conda deactivate
conda activate clair3


REF=/home/vbarbo/project_2021/paper_analysis/reference/genome/GRCh38.p13_all_chr.fasta
INPUT_BAM_DIR=/home/vbarbo/project_2021/paper_analysis/secondary_analyses/nanopore/files/ENCFF961HLO//data_manipulation
OUTPUT_DIR=/home/vbarbo/project_2021/paper_analysis/secondary_analyses/nanopore/files/ENCFF961HLO//variant_calling_from_nanopore/clair3
THREADS=10
PATH_TO_REPO=/home/vbarbo/project_2021/projects/lrRNAseqVariantCalling/tools/



### Guppy3,4-hac model downloaded from http://www.bio8.cs.hku.hk/clair3/clair3_models/r941_prom_hac_g360+g422.tar.gz
### See https://github.com/HKU-BAL/Clair3#pre-trained-models



### Clair3 (alone)
# call variants
/home/vbarbo/programs/Clair3/run_clair3.sh \
  --bam_fn=${INPUT_BAM_DIR}/aln_s.bam \
  --ref_fn=${REF} \
  --threads=${THREADS} \
  --platform="ont" \
  --model_path="/home/vbarbo/programs/Clair3/models/ont_guppy3-4" \
  --output=$OUTPUT_DIR/c3
## pileup
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
## full_alignment
# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/c3/full_alignment.vcf.gz \
  > $OUTPUT_DIR/c3/full_alignment_pass.vcf
# compress and index
bgzip -c $OUTPUT_DIR/c3/full_alignment_pass.vcf \
  > $OUTPUT_DIR/c3/full_alignment_pass.vcf.gz
bcftools index $OUTPUT_DIR/c3/full_alignment_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/c3/full_alignment_pass.vcf.gz
rm $OUTPUT_DIR/c3/full_alignment_pass.vcf






### SNCR + Clair3
/home/vbarbo/programs/Clair3/run_clair3.sh \
  --bam_fn=${INPUT_BAM_DIR}/aln_sncr.bam \
  --ref_fn=${REF} \
  --threads=${THREADS} \
  --platform="ont" \
  --model_path="/home/vbarbo/programs/Clair3/models/ont_guppy3-4" \
  --output=${OUTPUT_DIR}/c3_sncr
## pileup
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
## full_alignment
# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/c3_sncr/full_alignment.vcf.gz \
  > $OUTPUT_DIR/c3_sncr/full_alignment_pass.vcf
# compress and index
bgzip -c $OUTPUT_DIR/c3_sncr/full_alignment_pass.vcf \
  > $OUTPUT_DIR/c3_sncr/full_alignment_pass.vcf.gz
bcftools index $OUTPUT_DIR/c3_sncr/full_alignment_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/c3_sncr/full_alignment_pass.vcf.gz
rm $OUTPUT_DIR/c3_sncr/full_alignment_pass.vcf



### SNCR + FC + Clair3
/home/vbarbo/programs/Clair3/run_clair3.sh \
  --bam_fn=${INPUT_BAM_DIR}/aln_sncr_fc.bam \
  --ref_fn=${REF} \
  --threads=${THREADS} \
  --platform="ont" \
  --model_path="/home/vbarbo/programs/Clair3/models/ont_guppy3-4" \
  --output=${OUTPUT_DIR}/c3_sncr_fc
## pileup
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
## full_alignment
# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/c3_sncr_fc/full_alignment.vcf.gz \
  > $OUTPUT_DIR/c3_sncr_fc/full_alignment_pass.vcf
# compress and index
bgzip -c $OUTPUT_DIR/c3_sncr_fc/full_alignment_pass.vcf \
  > $OUTPUT_DIR/c3_sncr_fc/full_alignment_pass.vcf.gz
bcftools index $OUTPUT_DIR/c3_sncr_fc/full_alignment_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/c3_sncr_fc/full_alignment_pass.vcf.gz
rm $OUTPUT_DIR/c3_sncr_fc/full_alignment_pass.vcf




### delete temporary files
rm -r \
  $OUTPUT_DIR/c3/log/ \
  $OUTPUT_DIR/c3/tmp/ \
  $OUTPUT_DIR/c3_sncr/log/ \
  $OUTPUT_DIR/c3_sncr/tmp/ \
  $OUTPUT_DIR/c3_sncr_fc/log/ \
  $OUTPUT_DIR/c3_sncr_fc/tmp/

