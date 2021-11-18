#!/bin/bash


# use wtc-11 Iso-Seq data, gloria's data. downloaded from:
# https://sheynkman-lab-lifebit.s3.eu-west-1.amazonaws.com/Sample-Data/WTC11/XGSUV_20200804_S64049_PL100158447-1_A01.ccs.fastq.gz
# the data are ccs in fastq format.


# wtc-11
# iso-seq data
# duplicates not marked/removed
# only primary alignments



### inputs
CCS_BAM=/home/vbarbo/project_2021/datasets/wtc11/data/XGSUV_20200804_S64049_PL100158447-1_A01.ccs.bam
PRIMERS_FASTA=/home/vbarbo/project_2021/datasets/wtc11/data/NEB_primers.fasta
SAMPLE_NAME=XGSUV_20200804_S64049_PL100158447-1_A01
REF=/home/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta
OUTPUT_DIR=/home/vbarbo/project_2021/datasets/wtc11/manipulate_data
THREADS=30


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



### align reads to the genome of reference
minimap2 -ax splice \
  -t $THREADS \
  --secondary=no \
  $REF \
  $OUTPUT_DIR/${SAMPLE_NAME}_flnc.fastq \
  | samtools view -bSh -F 2308 - \
  > $OUTPUT_DIR/aln.bam
#  -uf -C5 \  ####### i should have used these parameters

# sort and index
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

#  index
samtools index \
  -@ $THREADS \
  $OUTPUT_DIR/aln_sncr.bam




### flagCorrection
### in R
times <- Sys.time()
variantCallingFromIsoSeq::flagCorrection(
  input_bam="/home/vbarbo/project_2021/datasets/wtc11/manipulate_data/aln_s.bam",
  input_sncr_bam="/home/vbarbo/project_2021/datasets/wtc11/manipulate_data/aln_sncr.bam",
  output_bam="/home/vbarbo/project_2021/datasets/wtc11/manipulate_data/aln_sncr_fc.bam",
  threads=30
)
(times <- Sys.time() - times)
# Time difference of 1.91619 hours

# index
samtools index -@ $THREADS \
  $OUTPUT_DIR/aln_sncr_fc.bam

