# this script is to create the table that stores the data needed to draw  the plot for the homopolymer analysis

# we kept only sites that:
#   * short-read variant density: x <= 3                          (ignore positions)
#   * ignore intronic regions (isoSeq_coverage > 0)               (ignore positions)
#   * DeepVariant QUAL: x >= 15                                   (remove method calling)
#   * short-read coverare : quantil 5% <= x <= 95%                (ignore positions)

#   * iso-seq coverage is higher than 20 reads;
#   * iso-seq variant density: x <= 3 (consider calls from all methods/pipelines);
#   * are far from splice junction (further than 20 bp);
#   * are classified as heterozygous alternative
#   * indels


### load packages
library(lrRNAseqBenchmark)
library(dplyr)
library(Biostrings)
library(tibble)
library(tidyr)
library(ggplot2)


### filter master table
# METHOD_NAMES <- c("dv_s_fc", "c3_mix", "gatk_s")
TRUTH_NAME <- "allen"
WTC11_MASTER_TABLE <- "/home/vbarbo/project_2021/paper_analysis/extra_files/master_tables/mt_wtc11_allMethods_filtered_v7_filterVariantDensityIsoSeqAllMethods.RData"
HOMOPOLYMERS_IN_REF <- "/home/vbarbo/project_2021/paper_analysis/extra_files/homopolymers_GRCh38.p13_all_chr.RData"
REF_FASTA <- "/home/vbarbo/project_2021/paper_analysis/reference/genome/GRCh38.p13_all_chr.fasta"
TRUTH_VCF_FILE <- "/home/vbarbo/project_2021/paper_analysis/wtc11/ground_truth/3546dc62_AH77TTBBXX_DS-229105_GCCAAT_recalibrated_subsetChromosomes_pass.vcf.gz"
METHOD_VCF_FILES <- c(
  "/home/vbarbo/project_2021/paper_analysis/wtc11/variant_calling_from_isoseq/deepvariant/dv_sncr_fc/deepvariant_calls_pass.vcf.gz",
  "/home/vbarbo/project_2021/paper_analysis/wtc11/variant_calling_from_isoseq/clair3/mix/pileup_pass_mix_nodup.recode.vcf.gz",
  "/home/vbarbo/project_2021/paper_analysis/wtc11/variant_calling_from_isoseq/gatk/isoSeq_wtc11.recal_pass.vcf.gz"
)


### get master table and filter it
k <- load(WTC11_MASTER_TABLE)
dat1 <- get(k)

dat1 <- filter(dat1, is_near_ss==0)


### get indels in homopolymers (homopolymer length equal to 1 means non homopolymer)
### also, filter sites by minimum iso-seq read coverage
load(HOMOPOLYMERS_IN_REF)
ref_fasta_seqs <- readDNAStringSet(REF_FASTA)
names(ref_fasta_seqs) <- sub(" .+", "", names(ref_fasta_seqs))

# dv_s_fc
dat_hom_dvSFc <- method_homopolymer_indels(input_table=dat1,
                                           first_method_name=TRUTH_NAME,
                                           second_method_name="dv_s_fc",
                                           vcf_first=TRUTH_VCF_FILE,
                                           vcf_second=METHOD_VCF_FILES[1],
                                           homopolymers=homopolymers,
                                           ref_fasta_seqs=ref_fasta_seqs,
                                           min_isoseq_coverage=20,
                                           genotyped_alt="find")
# c3_mix
dat_hom_c3Mix <- method_homopolymer_indels(input_table=dat1,
                                           first_method_name=TRUTH_NAME,
                                           second_method_name="c3_mix",
                                           vcf_first=TRUTH_VCF_FILE,
                                           vcf_second=METHOD_VCF_FILES[2],
                                           homopolymers=homopolymers,
                                           ref_fasta_seqs=ref_fasta_seqs,
                                           min_isoseq_coverage=20,
                                           genotyped_alt="find")
# gatk_s
dat_hom_gatkS <- method_homopolymer_indels(input_table=dat1,
                                           first_method_name=TRUTH_NAME,
                                           second_method_name="gatk_s",
                                           vcf_first=TRUTH_VCF_FILE,
                                           vcf_second=METHOD_VCF_FILES[3],
                                           homopolymers=homopolymers,
                                           ref_fasta_seqs=ref_fasta_seqs,
                                           min_isoseq_coverage=20,
                                           genotyped_alt="same")



### for each method to compare, extract the homopolymer information
### and concatenate the tables generating a single list
dat_full <-list(class_counts=NULL, dat_text=NULL)

# # # # # # # # # # # # # #
# # sncr+fc+deepvariant # #
# # # # # # # # # # # # # #

# deletion
k <- make_homopolymer_table_to_plot(
  input_hom_table <- dat_hom_dvSFc,
  variant_type <- "deletion",
  method_name <- "dv_s_fc",
  truth_name <- "allen",
  hom_length_intervals <- c(1,2,3,5,11,16,21),
  interval_names <- c("non-hp",2,"3-4","5-10", "11-15", "16-20", ">=21"),
  to_calculate <- "pre_rec_f1",
  output_method_name <- "SNCR+FC+DeepVariant"
)
dat_full <- mapply(rbind, dat_full, k, SIMPLIFY=FALSE)

# insertion
k <- make_homopolymer_table_to_plot(
  input_hom_table <- dat_hom_dvSFc,
  variant_type <- "insertion",
  method_name <- "dv_s_fc",
  truth_name <- "allen",
  hom_length_intervals <- c(1,2,3,5,11,16,21),
  interval_names <- c("non-hp",2,"3-4","5-10", "11-15", "16-20", ">=21"),
  to_calculate <- "pre_rec_f1",
  output_method_name <- "SNCR+FC+DeepVariant"
)
dat_full <- mapply(rbind, dat_full, k, SIMPLIFY=FALSE)

# # # # # # # # # #
# # clair3 mix  # #
# # # # # # # # # #

k <- make_homopolymer_table_to_plot(
  input_hom_table <- dat_hom_c3Mix,
  variant_type <- "deletion",
  method_name <- "c3_mix",
  truth_name <- "allen",
  hom_length_intervals <- c(1,2,3,5,11,16,21),
  interval_names <- c("non-hp",2,"3-4","5-10", "11-15", "16-20", ">=21"),
  to_calculate <- "pre_rec_f1",
  output_method_name <- "Clair3 mix"
)
dat_full <- mapply(rbind, dat_full, k, SIMPLIFY=FALSE)

k <- make_homopolymer_table_to_plot(
  input_hom_table <- dat_hom_c3Mix,
  variant_type <- "insertion",
  method_name <- "c3_mix",
  truth_name <- "allen",
  hom_length_intervals <- c(1,2,3,5,11,16,21),
  interval_names <- c("non-hp",2,"3-4","5-10", "11-15", "16-20", ">=21"),
  to_calculate <- "pre_rec_f1",
  output_method_name <- "Clair3 mix"
)
dat_full <- mapply(rbind, dat_full, k, SIMPLIFY=FALSE)

# # # # # # # # #
# # sncr+gatk # #
# # # # # # # # #

k <- make_homopolymer_table_to_plot(
  input_hom_table <- dat_hom_gatkS,
  variant_type <- "deletion",
  method_name <- "gatk_s",
  truth_name <- "allen",
  hom_length_intervals <- c(1,2,3,5,11,16,21),
  interval_names <- c("non-hp",2,"3-4","5-10", "11-15", "16-20", ">=21"),
  to_calculate <- "pre_rec_f1",
  output_method_name <- "SNCR+GATK"
)
dat_full <- mapply(rbind, dat_full, k, SIMPLIFY=FALSE)

k <- make_homopolymer_table_to_plot(
  input_hom_table <- dat_hom_gatkS,
  variant_type <- "insertion",
  method_name <- "gatk_s",
  truth_name <- "allen",
  hom_length_intervals <- c(1,2,3,5,11,16,21),
  interval_names <- c("non-hp",2,"3-4","5-10", "11-15", "16-20", ">=21"),
  to_calculate <- "pre_rec_f1",
  output_method_name <- "SNCR+GATK"
)
dat_full <- mapply(rbind, dat_full, k, SIMPLIFY=FALSE)



### some manipulations
k <- recode(levels(dat_full$class_counts$Classification),
            "precision"="Precision", "sensitivity"="Recall", "f1Score"="F1-score")
levels(dat_full$class_counts$Classification) <- k

dat_full$dat_text$y <- max(dat_full$dat_text$y)

k <- c("SNCR+FC+DeepVariant", "Clair3 mix", "SNCR+GATK")
dat_full$class_counts$method <- factor(dat_full$class_counts$method,
                                       levels=k, ordered=TRUE)
dat_full$dat_text$method <- factor(dat_full$dat_text$method,
                                   levels=k, ordered=TRUE)


dat_homopolymer <- dat_full

### save object `dat_homopolymer` to a file.
### this is the data used to draw the chart for the homopolymer analysis
save(dat_homopolymer, file="/home/vbarbo/project_2021/paper_analysis/extra_files/dat_homopolymer_analysis.RData")


