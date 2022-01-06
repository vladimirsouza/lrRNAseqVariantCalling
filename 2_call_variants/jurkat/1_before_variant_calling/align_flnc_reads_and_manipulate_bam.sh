#!/bin/bash



### inputs
FLNC_BAM=/home/vbarbo/project_2021/datasets/gloria_data/jurkat_data/jurkat.flnc.bam
REF=/home/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta
OUTPUT_ROOT_DIR=/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls
OUTPUT_DIR=/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/noMarkDuplicate
THREADS=30



### convert the unaligned bam file to fastq file
#flnc_fastq=${FLNC_BAM%bam}fastq
flnc_fastq=${OUTPUT_ROOT_DIR}/$(basename $FLNC_BAM)
flnc_fastq=${flnc_fastq%bam}fastq
bamToFastq \
  -i $FLNC_BAM \
  -fq $flnc_fastq



### align reads to the genome of reference
minimap2 -ax splice \
  -t $THREADS \
  --secondary=no \
  $REF \
  $flnc_fastq \
  | samtools view -bSh -F 2308 - \
  > $OUTPUT_ROOT_DIR/aln.bam



### sort and index
samtools sort \
  -@ $THREADS \
  -o $OUTPUT_ROOT_DIR/aln.bam \
  $OUTPUT_ROOT_DIR/aln.bam

samtools index \
  -@ $THREADS \
  $OUTPUT_ROOT_DIR/aln.bam







### Split N-cigar reads - remove intronic regions - each exon is a read
### there is no need to reassign the mapping quality:
### (https://lh3.github.io/minimap2/minimap2.html):
### "Mapping quality (0-255 with 255 for missing)"
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" SplitNCigarReads \
  -R $REF \
  -I $OUTPUT_ROOT_DIR/aln.bam \
  -O $OUTPUT_DIR/aln_split.bam



### flagCorrection
### in R
times <- Sys.time()
library(variantCallingFromIsoSeq)
flagCorrection(
  input_bam="/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/aln.bam",
  input_sncr_bam="/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/noMarkDuplicate/aln_split.bam",
  output_bam="/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/noMarkDuplicate/aln_split_flagCorrection.bam",
  threads=30
)
(times <- Sys.time() - times)
#Time difference of ??? mins



### index
samtools index -@ $THREADS $OUTPUT_DIR/aln_split_flagCorrection.bam


