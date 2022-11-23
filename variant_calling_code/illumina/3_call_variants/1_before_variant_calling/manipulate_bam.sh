#!/bin/bash



# this script generates all different BAM files (using or not SNCR/flagCorrection) to call variants from Illumina RNA-Seq short reads (WTC-11 dataset).
# we use the same data used for the ASE analysis.

# The imput BAM was generated in script '1_count_reads_per_allele.sh'. Non-primoary alignments (unmapped, secondary, supplementary) were removed.
# Duplicates were kept. Read groups were added.




### inputs
SHORT_READ_BAM=/home/vbarbo/project_2021/paper_analysis/wtc11/ase_analysis/star_aln/GSE175048_Aligned.sortedByCoord.out_primary_readGroupAddedSameOrder.bam
OUTPUT_DIR=/home/vbarbo/project_2021/paper_analysis/secondary_analyses/short_reads/files/data_manipulation
REF=/home/vbarbo/project_2021/paper_analysis/reference/genome/GRCh38.p13_all_chr.fasta
THREADS=30
PATH_TO_FC=/home/vbarbo/project_2021/projects/lrRNAseqVariantCalling/tools/flagCorrection.r





### remove introns with SplitNCigarReads
sncr_bam_base_name=$( basename ${SHORT_READ_BAM%.bam}_sncr.bam )
mkdir -p ${OUTPUT_DIR}/temp
gatk --java-options "-Xmx10G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" SplitNCigarReads \
    -R ${REF} \
    -I ${SHORT_READ_BAM} \
    -O ${OUTPUT_DIR}/${sncr_bam_base_name} \
    --tmp-dir ${OUTPUT_DIR}/temp



###  index
samtools index \
  -@ ${THREADS} \
  ${OUTPUT_DIR}/${sncr_bam_base_name}



### flagCorrection
fc_bam_base_name=${sncr_bam_base_name%.bam}_fc.bam
Rscript ${PATH_TO_FC} \
  ${SHORT_READ_BAM} \
  ${OUTPUT_DIR}/${sncr_bam_base_name} \
  ${OUTPUT_DIR}/${fc_bam_base_name} \
  ${THREADS}



### index
samtools index \
  -@ ${THREADS} \
  ${OUTPUT_DIR}/${fc_bam_base_name}


