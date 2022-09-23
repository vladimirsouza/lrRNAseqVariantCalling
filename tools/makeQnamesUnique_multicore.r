#!/usr/bin/env Rscript


# When GATK's SplitNCigarReads function splits exons of a read, it makes 
# the qnames of all new reads to the be same of the original one. This 
# script makes all qnames of a bam file to be unique, by adding '_<int>',
# at the end of the each qname, where <int> is a integer coresponding to
# the number of the exon. Exemplos:
# 
# read_number     exon_number    qname_of_input_bam     qname_of_output_bam
#           1               1           read1_qname           read1_qname_1
#           1               2           read1_qname           read1_qname_2
#           2               1           read2_qname           read2_qname_1
#           2               2           read2_qname           read2_qname_2
# 
#
# The input parameters are (following the right order):
# * `FAULTY_BAM`, path to input BAM file;
# * `OUTPUT_BAM`, path to output BAM file;
# * `THREADS`, number of cores to use.
#
#
# To run makeQnamesUnique:
# Rscript $PATH_TO_REPO/makeQnamesUnique.r \
#   $INPUT_BAM \
#   $OUTPUT_BAM \
#   $THREADS




### input variables
INPUT_BAM=commandArgs(TRUE)[1]
OUTPUT_BAM=commandArgs(TRUE)[2]
THREADS=as.integer(commandArgs(TRUE)[3])




### time when of start
time1 <- Sys.time()


### test inputs
if( !file.exists(INPUT_BAM) ) stop( gettextf("File %s doesn't exist", INPUT_BAM) )
output_dir <- dirname(OUTPUT_BAM)
if( !dir.exists(output_dir) ) stop( gettextf("Directory %s doesn't exist", output_dir) )


### libraries
suppressPackageStartupMessages( library(foreach) )
suppressPackageStartupMessages( library(doParallel) )


### create root temp dir
temp_dir <- tempfile(pattern="temp_dir_", tmpdir=output_dir)
dir.create(temp_dir)


### write sam header from bam
sam_header <- tempfile(pattern="sam_header_", tmpdir=temp_dir, fileext=".sam")
cmd <- gettextf("samtools view -H %s > %s", INPUT_BAM, sam_header)
system(cmd)


### get chromosome names
cmd <- gettextf("perl -lne '/SN:(\\S+)/ and print $1' %s", sam_header)
chr_names <- system(cmd, intern=TRUE)


### define some variables
# temp dir for each chromosome
chr_temp_dirs <- file.path(temp_dir, chr_names)
# faulty sam body files
faulty_sam_files <- tempfile(pattern="faulty_sam_file_", tmpdir=chr_temp_dirs, fileext=".sam")
# fixed qnames txt files
fixed_qname_files <- tempfile(pattern="fixed_qname_", tmpdir=chr_temp_dirs, fileext=".txt")
# fixed sam body files
fixed_sam_files <- tempfile(pattern="fixed_sam_files_", tmpdir=chr_temp_dirs, fileext=".sam")
# headered fixed sam files
headered_fixed_sam_files <- tempfile(pattern="headered_fixed_sam_files_", tmpdir=chr_temp_dirs, fileext=".sam")
# fixed bam files
bam_files <- gettextf("%s/%s.bam", temp_dir, chr_names)


### make fixed-qnames sam files for each chromosome
make_qnames_unique <- function(chr_temp_dirs_i, chr_names_i, faulty_sam_files_i,
                               fixed_qname_files_i, fixed_sam_files_i,
                               headered_fixed_sam_files_i, bam_files_i,
                               sam_header){
  
  ### create a temp dir for the chromosome
  dir.create(chr_temp_dirs_i)
  
  ### create a SAM file for each chromosome
  cmd <- gettextf("samtools view %s %s > %s", INPUT_BAM, chr_names_i, faulty_sam_files_i)
  system(cmd)
  
  ### load qname
  cmd <- gettextf("cut -f 1 %s", faulty_sam_files_i)
  qnames <- system(cmd, intern=TRUE)
  
  ### fix qnames
  qnames_uni <- unique(qnames)
  invisible(
    lapply(qnames_uni, function(qnames_uni_i){
      # qnames_uni_i <- qnames_uni[[1]]
      
      k <- qnames == qnames_uni_i
      qnames[k] <<- paste0(qnames_uni_i, "_", 1:sum(k))
    })
  )
  
  ### save fixed qnames to a file
  write.table(qnames, file=fixed_qname_files_i, quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  ### replace faulty qnames with fixed ones
  cmd <- gettextf("awk 'FNR==NR{a[NR]=$1;next}{$1=a[FNR]}1' OFS='\t' %s %s > %s",
                  fixed_qname_files_i, faulty_sam_files_i, fixed_sam_files_i)
  system(cmd)
  
  ### add header to sam
  cmd <- gettextf("cat %s %s > %s",
                  sam_header, fixed_sam_files_i, headered_fixed_sam_files_i)
  system(cmd)
  
  ### convert sam to bam
  cmd <- gettextf("samtools view -S -b %s > %s",
                  headered_fixed_sam_files_i, bam_files_i)
  system(cmd)
  
  ### delete chromosome temp files and mark dir as done
  unlink(chr_temp_dirs_i, recursive=TRUE)
}


### maximun number of cores is the number of chomosomes
threads <- min( THREADS, length(chr_names) )


### run function using multiple cores
registerDoParallel(cores=threads)
invisible(
  foreach(chr_temp_dirs_i=chr_temp_dirs,
          chr_names_i=chr_names,
          faulty_sam_files_i=faulty_sam_files,
          fixed_qname_files_i=fixed_qname_files,
          fixed_sam_files_i=fixed_sam_files,
          headered_fixed_sam_files_i=headered_fixed_sam_files,
          bam_files_i=bam_files) %dopar%
    make_qnames_unique(chr_temp_dirs_i,
                       chr_names_i,
                       faulty_sam_files_i,
                       fixed_qname_files_i,
                       fixed_sam_files_i,
                       headered_fixed_sam_files_i,
                       bam_files_i,
                       sam_header)
)


### if `OUTPUT_BAM` doesn't contain the bam extension, add it
if( ! grepl("\\.bam$", OUTPUT_BAM) )
  OUTPUT_BAM <- paste0(OUTPUT_BAM, ".bam")


### merge bams
all_bams_to_merge <- paste(bam_files, collapse=" ")
cmd <- gettextf("samtools merge -c -p -f -@ %i %s %s", THREADS, OUTPUT_BAM, all_bams_to_merge)
system(cmd)


### delete temp dir
unlink(temp_dir, recursive=TRUE)


### report time spent
time2 <- Sys.time()
Time <- time2 - time1
Time <- round( as.numeric(Time, units = "mins"), 2 )
cat( gettextf("\nmakeQnamesUnique finished after %s minutes.\n\n", Time) )

