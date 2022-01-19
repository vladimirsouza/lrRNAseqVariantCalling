#!/usr/bin/env Rscript


# This script takes the FLAG of each read in `ORI_BAM` and apply them to reads of the same QNAME in `SNCR_BAM`.
#
# The input parameters are (following the right order):
# * `ORI_BAM`, path of the original BAM file, BAM before spliting reads with GATK's SplitNCigarReads (SNCR) function;
# * `SNCR_BAM`, path of the split BAM file, BAM after spliting reads with SNCR;
# * `OUTPUT_DIR`, path of the directory to write the output BAM file;
# * `OUTPUT_BAM_NAME`, name of the output BAM file, the '.bam' extension is added in case it is missing; and
# * `THREADS`, number of cores to use.
#
# To run flagCorrection:
# Rscript $PATH_TO_FLAGCORRECTION/flagCorrection.r \
#   $ORI_BAM \
#   $SNCR_BAM \
#   $OUTPUT_DIR \
#   $OUTPUT_BAM_NAME \
#   $THREADS


time1 <- Sys.time()


ORI_BAM=commandArgs(TRUE)[1]
SNCR_BAM=commandArgs(TRUE)[2]
OUTPUT_DIR=commandArgs(TRUE)[3]
OUTPUT_BAM_NAME=commandArgs(TRUE)[4]
THREADS=as.integer(commandArgs(TRUE)[5])


### load packages
library(Rsamtools)
library(foreach)
library(doParallel)


### test inputs
if( !file.exists(ORI_BAM) ) stop( gettextf("File %s doesn't exist", ORI_BAM) )
if( !file.exists(SNCR_BAM) ) stop( gettextf("File %s doesn't exist", SNCR_BAM) )
if( !dir.exists(OUTPUT_DIR) ) stop( gettextf("Directory %s doesn't exist", OUTPUT_DIR) )


### create temp dir
temp_dir <- tempfile("flagCorrection_temp_dir_", tmpdir=OUTPUT_DIR)
dir.create(temp_dir)


### get flags and qnames from the original bam
k <- scanBam( ORI_BAM, param=ScanBamParam(what=c("qname", "flag")) ) [[1]]
ori_bam_qname_flag <- setNames(k$flag, k$qname)


### get headers only
bam_header <- file.path(temp_dir, "header.sam")
cmd <- gettextf("samtools view -H %s > %s", ORI_BAM, bam_header)
system(cmd)


### define some variables
# chromosome names
cmd <- gettextf("samtools view -H %s | perl -lne '/SN:(\\S+)/ and print $1'", SNCR_BAM)
chr_names <- system(cmd, intern=TRUE)
# split sam addresses
k <- file.path(temp_dir, chr_names)
split_sncr_sams <- paste0(k, ".sam")
# corrected flags
corrected_flags_files <- paste0(k, "_corrected_flags.txt")
# corrected split sam addresses
corrected_split_sncr_sams <- paste0(k, "_corrected.sam")
# headered corrected split sam addresses
headered_corrected_split_sncr_sams <- paste0(k, "_corrected_headered.sam")
# headered corrected split bam addresses
final_bams <- paste0(k, "_final.bam")
# threads
threads <- min( THREADS, length(chr_names) )


### function to split bam by chromosomes
correct_flags_per_chromosome <- function(chr_names_i, split_sncr_sams_i, corrected_flags_files_i, corrected_split_sncr_sams_i,
                                         headered_corrected_split_sncr_sams_i, final_bams_i, SNCR_BAM, ori_bam_qname_flag, bam_header){
  
  # create the SAM file of the chromosome
  cmd <- gettextf("samtools view %s %s > %s", SNCR_BAM, chr_names_i, split_sncr_sams_i)
  system(cmd)
  
  # load QNAMEs from sncr bam
  cmd <- gettextf("cut -f 1 %s", split_sncr_sams_i)
  qname_sncr <- system(cmd, intern=TRUE)
  
  # get flags corrected
  corrected_flag <- unname( ori_bam_qname_flag[qname_sncr] )
  write.table(corrected_flag, file=corrected_flags_files_i, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  
  # replace the flags of sam
  cmd <- gettextf("awk 'FNR==NR{a[NR]=$1;next}{$2=a[FNR]}1' OFS='\t' %s %s > %s", corrected_flags_files_i, split_sncr_sams_i, corrected_split_sncr_sams_i)
  system(cmd)
  unlink(split_sncr_sams_i)
  
  # add headers to sam
  cmd <- gettextf("cat %s %s > %s", bam_header, corrected_split_sncr_sams_i, headered_corrected_split_sncr_sams_i)
  system(cmd)
  unlink(corrected_split_sncr_sams_i)
  
  # sam to bam
  cmd <- gettextf("samtools view -S -b %s > %s", headered_corrected_split_sncr_sams_i, final_bams_i)
  system(cmd)
  unlink(headered_corrected_split_sncr_sams_i)
}


### get corrected SAMs
registerDoParallel(cores=threads)
result <- foreach(chr_names_i = chr_names,
                  split_sncr_sams_i = split_sncr_sams,
                  corrected_flags_files_i = corrected_flags_files,
                  corrected_split_sncr_sams_i = corrected_split_sncr_sams,
                  headered_corrected_split_sncr_sams_i = headered_corrected_split_sncr_sams,
                  final_bams_i = final_bams) %dopar%
  correct_flags_per_chromosome(chr_names_i,
                               split_sncr_sams_i,
                               corrected_flags_files_i,
                               corrected_split_sncr_sams_i,
                               headered_corrected_split_sncr_sams_i,
                               final_bams_i,
                               SNCR_BAM,
                               ori_bam_qname_flag,
                               bam_header)
unlink( c(corrected_flags_files, bam_header) )


### output bam path
output_bam <- file.path(OUTPUT_DIR, OUTPUT_BAM_NAME)
if( !grepl("\\.bam$", OUTPUT_BAM_NAME) ){
  output_bam <- paste0(output_bam, ".bam")
}


### merge bams
all_bams_to_merge <- paste(final_bams, collapse=" ")
cmd <- gettextf("samtools merge -c -p -f -@ %i %s %s", THREADS, output_bam, all_bams_to_merge)
system(cmd)


### delete temp files/dir
unlink(temp_dir, recursive=TRUE)


### report time spent
time2 <- Sys.time()
Time <- time2 - time1
Time <- round( as.numeric(Time, units = "mins"), 2 )
cat( gettextf("\nflagCorrection finished after %s minutes.\n\n", Time) )


