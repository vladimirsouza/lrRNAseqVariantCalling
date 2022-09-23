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
# * `OUTPUT_BAM`, path to output BAM file.
#
#
# To run makeQnamesUnique:
# Rscript $PATH_TO_REPO/makeQnamesUnique.r \
#   $INPUT_BAM \
#   $OUTPUT_BAM





### input variables
INPUT_BAM=commandArgs(TRUE)[1]
OUTPUT_BAM=commandArgs(TRUE)[2]



### test inputs
if( !file.exists(INPUT_BAM) ) stop( gettextf("File %s doesn't exist", INPUT_BAM) )
output_dir <- dirname(OUTPUT_BAM)
if( !dir.exists(output_dir) ) stop( gettextf("Directory %s doesn't exist", output_dir) )


### libraries
suppressMessages( library(Rsamtools) )


### load qnames
qn <- scanBam(INPUT_BAM, param=ScanBamParam(what="qname")) [[1]] $qname


### fix qnames
qn_uni <- unique(qn)
invisible(
  lapply(qn_uni, function(qn_uni_i){
    # qn_uni_i <- qn_uni[[1]]
    
    k <- qn == qn_uni_i
    qn_new <- paste0(qn_uni_i, "_", 1:sum(k))
    qn[k] <<- qn_new
  })
)


### save the fixed qnames to a file
temp_dir <- tempfile(pattern="temp_dir_", tmpdir=output_dir)
dir.create(temp_dir)
fixed_qnames_file <- tempfile(pattern="fixed_qnames_", tmpdir=temp_dir, fileext=".txt")
write.table(qn, file=fixed_qnames_file, quote=FALSE, row.names=FALSE, col.names=FALSE)


### write sam header from bam
faulty_sam_header <- tempfile(pattern="faulty_sam_header_", tmpdir=temp_dir, fileext=".sam")
cmd <- gettextf("samtools view -H %s > %s", INPUT_BAM, faulty_sam_header)
system(cmd)


### write sam body from bam
faulty_sam_body <- tempfile(pattern="faulty_sam_body_", tmpdir=temp_dir, fileext=".sam")
cmd <- gettextf("samtools view %s > %s", INPUT_BAM, faulty_sam_body)
system(cmd)


### replace the faulty qnames with the fixed ones
fixed_sam_body <- tempfile(pattern="fixed_sam_body_", tmpdir=temp_dir, fileext=".sam")
cmd <- gettextf("awk 'FNR==NR{a[NR]=$1;next}{$1=a[FNR]}1' OFS='\t' %s %s > %s",
                fixed_qnames_file, faulty_sam_body, fixed_sam_body)
system(cmd)


### join header and body
fixed_sam <- tempfile(pattern="fixed_sam_", tmpdir=temp_dir, fileext=".sam")
cmd <- gettextf("cat %s %s > %s", faulty_sam_header, fixed_sam_body, fixed_sam)
system(cmd)


### if `OUTPUT_BAM` doesn't contain the bam extension, add it
if( ! grepl("\\.bam$", OUTPUT_BAM) )
  OUTPUT_BAM <- paste0(OUTPUT_BAM, ".bam")


### sam to the final bam file \o/
cmd <- gettextf("samtools view -h -b %s > %s", fixed_sam, OUTPUT_BAM)
system(cmd)


### delete temp files/dir
unlink(temp_dir, recursive=TRUE)


