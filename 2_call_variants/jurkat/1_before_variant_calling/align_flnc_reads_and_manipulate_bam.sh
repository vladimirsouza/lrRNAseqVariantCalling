#!/bin/bash



### inputs
FLNC_BAM=/home/vbarbo/project_2021/paper_analysis/jurkat/data/iso_seq/jurkat.flnc.bam
REF=/home/vbarbo/project_2021/paper_analysis/reference/genome/GRCh38.p13_all_chr.fasta
OUTPUT_DIR=/home/vbarbo/project_2021/paper_analysis/jurkat/data_manipulation
THREADS=30
PATH_TO_FC=/home/vbarbo/project_2021/projects/lrRNA-seq_variant_calling/flagCorrection.r



### convert the unaligned bam file to fastq file
bamToFastq \
  -i ${FLNC_BAM} \
  -fq ${OUTPUT_DIR}/jurkat_flnc.fastq



### align reads to the genome of reference
minimap2 -ax splice \
  -uf -C5 \
  -t $THREADS \
  --secondary=no \
  $REF \
  ${OUTPUT_DIR}/jurkat_flnc.fastq \
  | samtools view -bSh -F 2308 - \
  > $OUTPUT_DIR/aln.bam



### sort and index
samtools sort \
  -@ $THREADS \
  -o $OUTPUT_DIR/aln.bam \
  $OUTPUT_DIR/aln.bam

samtools index \
  -@ $THREADS \
  $OUTPUT_DIR/aln.bam




### Split N-cigar reads - remove intronic regions - each exon is a read
### there is no need to reassign the mapping quality:
### (https://lh3.github.io/minimap2/minimap2.html):
### "Mapping quality (0-255 with 255 for missing)"
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" SplitNCigarReads \
  -R $REF \
  -I $OUTPUT_DIR/aln.bam \
  -O $OUTPUT_DIR/aln_sncr.bam



### flagCorrection
Rscript $PATH_TO_FC \
  $OUTPUT_DIR/aln.bam \
  $OUTPUT_DIR/aln_sncr.bam \
  $OUTPUT_DIR/aln_sncr_fc.bam \
  $THREADS


### index
samtools index \
  -@ $THREADS \
  $OUTPUT_DIR/aln_sncr_fc.bam


