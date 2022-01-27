#!/bin/bash



# this script generates all different BAM files (using or not SNCR/flagCorrection) to call variants from Iso-Seq data (WTC-11 dataset).



# wtc-11
# iso-seq data
# duplicates not marked/removed
# only primary alignments



### inputs
CCS_BAM=/home/vbarbo/project_2021/paper_analysis/wtc11/data/iso_seq/XGSUV_20200804_S64049_PL100158447-1_A01.ccs.bam
PRIMERS_FASTA=/home/vbarbo/project_2021/paper_analysis/wtc11/data/iso_seq/NEB_primers.fasta
SAMPLE_NAME=XGSUV_20200804_S64049_PL100158447-1_A01
REF=/home/vbarbo/project_2021/paper_analysis/reference/genome/GRCh38.p13_all_chr.fasta
OUTPUT_DIR=/home/vbarbo/project_2021/paper_analysis/wtc11/data_manipulation
THREADS=30
PATH_TO_FC=/home/vbarbo/project_2021/projects/lrRNA-seq_variant_calling/flagCorrection.r


### isoseq pipeline: get flnc reads
conda deactivate
conda activate isoSeq3


# Remove primers
lima \
  -j $THREADS \
  $CCS_BAM \
  $PRIMERS_FASTA \
  $OUTPUT_DIR/${SAMPLE_NAME}_fl.bam \
  --isoseq


# refine
isoseq3 refine \
  $OUTPUT_DIR/${SAMPLE_NAME}_fl.NEB_5p--NEB_3p.bam \
  $PRIMERS_FASTA \
  $OUTPUT_DIR/${SAMPLE_NAME}_flnc.bam \
  --require-polya


### convert the unaligned bam file (flnc) to fastq file
conda deactivate
conda activate base

bamToFastq \
  -i $OUTPUT_DIR/${SAMPLE_NAME}_flnc.bam \
  -fq $OUTPUT_DIR/${SAMPLE_NAME}_flnc.fastq


### align reads to the genome of reference and remove secondary and supplementary alignments (keep duplicates)
minimap2 -ax splice \
  -uf -C5 \
  -t $THREADS \
  --secondary=no \
  $REF \
  $OUTPUT_DIR/${SAMPLE_NAME}_flnc.fastq \
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


### Split N-cigar reads, remove intronic regions, make one read for each exon
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" SplitNCigarReads \
  -R $REF \
  -I $OUTPUT_DIR/aln_s.bam \
  -O $OUTPUT_DIR/aln_sncr.bam


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


