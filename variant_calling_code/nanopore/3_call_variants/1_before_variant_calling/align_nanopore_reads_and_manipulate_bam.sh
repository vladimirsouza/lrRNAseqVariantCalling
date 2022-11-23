#!/bin/bash



# this script generates all different BAM files (using or not SNCR/flagCorrection) to call variants from Nanopore RNA-Seq raw sequencing data (WTC-11 dataset).



# wtc-11
# Nanopore RNA-Seq raw sequencing data
# duplicates not marked/removed
# only primary alignments



### inputs
FASTQ=/home/vbarbo/project_2021/paper_analysis/secondary_analyses/nanopore/datasets/ENCSR539ZXJ/ENCFF961HLO.fastq.gz
SAMPLE_NAME=ENCFF961HLO
REF=/home/vbarbo/project_2021/paper_analysis/reference/genome/GRCh38.p13_all_chr.fasta
OUTPUT_DIR=/home/vbarbo/project_2021/paper_analysis/secondary_analyses/nanopore/files/ENCFF961HLO/data_manipulation
THREADS=30
PATH_TO_FC=/home/vbarbo/project_2021/projects/lrRNAseqVariantCalling/tools/flagCorrection.r




### align reads to the genome of reference and remove secondary and supplementary alignments (keep duplicates)
minimap2 -ax splice \
  -uf -k14 \
  -t $THREADS \
  --secondary=no \
  $REF \
  $FASTQ \
  | samtools view -bSh -F 2308 - \
  > $OUTPUT_DIR/aln.bam



### sort and index
samtools sort \
  -@ $THREADS \
  -o $OUTPUT_DIR/aln_s.bam \
  $OUTPUT_DIR/aln.bam
samtools index \
  -@ $THREADS \
  $OUTPUT_DIR/aln_s.bam



### to avaid java complaning of heap space, split the BAM file by chromosome and use SplitNCigarReads separately on each BAM

# split BAM by chromosome
mkdir -p ${OUTPUT_DIR}/split_bams/temp
bamtools split \
  -in $OUTPUT_DIR/aln_s.bam \
  -stub ${OUTPUT_DIR}/split_bams/aln_s \
  -reference

# split N-cigar reads
split_bam_s=( $(find ${OUTPUT_DIR}/split_bams/aln_s*) )
for i in ${split_bam_s[@]}
do
  gatk --java-options "-Xmx10G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" SplitNCigarReads \
    -R $REF \
    -I ${i} \
    -O ${i%.bam}_sncr.bam \
    --tmp-dir ${OUTPUT_DIR}/split_bams/temp
done
rmdir ${OUTPUT_DIR}/split_bams/temp

# merge BAMs
all_chr_bams=$(find ${OUTPUT_DIR}/split_bams -name "*_sncr.bam" )
samtools merge \
  ${OUTPUT_DIR}/aln_sncr.bam \
  ${all_chr_bams}
rm -r ${OUTPUT_DIR}/split_bams



###  index
samtools index \
  -@ $THREADS \
  $OUTPUT_DIR/aln_sncr.bam



### flagCorrection
Rscript $PATH_TO_FC \
  $OUTPUT_DIR/aln_s.bam \
  $OUTPUT_DIR/aln_sncr.bam \
  $OUTPUT_DIR/aln_sncr_fc.bam \
  $THREADS



### index
samtools index \
  -@ $THREADS \
  $OUTPUT_DIR/aln_sncr_fc.bam


