# this script is to create the table that stores the data needed to draw plots for the homopolymer analysis

# besides the comomn filters, i filterd out variants that:
#   * iso-seq coverage is lower than 20 reads;
#   * are in regions (201 bp window) that the variant density, from iso-seq data, is >3 (consider calls from all methods/pipelines);
#   * are near to splice junction (not further than 20 bp);
#   * are classified as heterozygous alternative



### load packages
library(variantCallingFromIsoSeq)
library(dplyr)
library(Biostrings)
library(tibble)
library(tidyr)
library(ggplot2)


### filter master table
METHOD_NAMES <- c("dv", "dv_s", "dv_s_fc",
                  "c3", "c3_s", "c3_s_fc", "c3_mix",
                  "gatk_s")
TRUTH_NAME <- "allen"
load("~/project_2021/scripts/paper_1/master_tables/wtc11/mt_wtc11_allMethods_notFiltered_v5.RData")
dat1 <- mt_wtc11_allMethods_notFiltered

shortread_cover_quantiles <- quantile(dat1$shortRead_coverage, probs=c(.05, .95))

k <- !( dat1$variantDensity_allen > 3 )
k1 <- !( dat1$variantDensity_dv_dvS_dvSFc_c3_c3S_c3SFc_c3Mix_gatkS > 3 )
k <- k & k1
dat1 <- dat1[k,]

dat1 <- filter(dat1, isoSeq_coverage>0)
dim(dat1)

dat1 <- filter(dat1, shortRead_coverage >= shortread_cover_quantiles[1] &
                 shortRead_coverage <= shortread_cover_quantiles[2])

k <- which( dat1$qual_dv < 15 )
dat1$in_dv[k] <- 0
dat1 <- add_method_vs_truth_comparison_to_master_table(
  dat1,
  METHOD_NAMES[1],
  TRUTH_NAME,
  replace_column=TRUE
)

k <- which( dat1$qual_dvS < 15 )
dat1$in_dv_s[k] <- 0
dat1 <- add_method_vs_truth_comparison_to_master_table(
  dat1,
  METHOD_NAMES[2],
  TRUTH_NAME,
  replace_column=TRUE
)

k <- which( dat1$qual_dvSFc < 15 )
dat1$in_dv_s_fc[k] <- 0
dat1 <- add_method_vs_truth_comparison_to_master_table(
  dat1,
  METHOD_NAMES[3],
  TRUTH_NAME,
  replace_column=TRUE
)

dat1 <- filter(dat1, is_near_ss==0)




### get indels in homopolymers (homopolymer length equal to 1 means non homopolymer)
### also, filter sites by minimum iso-seq read coverage
load("~/load_later/homopolymers_GRCh38.p13.RData")
ref_fasta_seqs <- readDNAStringSet("/home/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta")
names(ref_fasta_seqs) <- sub(" .+", "", names(ref_fasta_seqs))
TRUTH_VCF_FILE <- "/home/vbarbo/project_2021/datasets/wtc11/truth/allen_institute_wtc11_shortreads_wgs/3546dc62_AH77TTBBXX_DS-229105_GCCAAT_recalibrated_subsetChromosomes_pass.vcf.gz"
METHOD_VCF_FILES <- c(
  "~/project_2021/datasets/wtc11/methods_to_comp/deepvariant/deepvariant_calls_pass.vcf.gz",
  "~/project_2021/datasets/wtc11/methods_to_comp/clair3/mix/pileup_pass_mix_nodup.recode.vcf.gz",
  "~/project_2021/datasets/wtc11/methods_to_comp/gatk/isoSeq_wtc11.recal_pass.vcf.gz"
)
# dv_s_fc
dat_hom_dvSFc <- method_homopolymer_indels(input_table=dat1,
                                           first_method_name="allen",
                                           second_method_name="dv_s_fc",
                                           vcf_first=TRUTH_VCF_FILE,
                                           vcf_second=METHOD_VCF_FILES[[1]],
                                           homopolymers=homopolymers,
                                           ref_fasta_seqs=ref_fasta_seqs,
                                           min_isoseq_coverage=20,
                                           genotyped_alt="find")
# c3_mix
dat_hom_c3Mix <- method_homopolymer_indels(input_table=dat1,
                                           first_method_name="allen",
                                           second_method_name="c3_mix",
                                           vcf_first=TRUTH_VCF_FILE,
                                           vcf_second=METHOD_VCF_FILES[[2]],
                                           homopolymers=homopolymers,
                                           ref_fasta_seqs=ref_fasta_seqs,
                                           min_isoseq_coverage=20,
                                           genotyped_alt="find")
# gatk_s
dat_hom_gatkS <- method_homopolymer_indels(input_table=dat1,
                                           first_method_name="allen",
                                           second_method_name="gatk_s",
                                           vcf_first=TRUTH_VCF_FILE,
                                           vcf_second=METHOD_VCF_FILES[[3]],
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


### save object `dat_full` to a file.
### this is the data used to draw the chart for the homopolymer analysis
save(dat_full, file="~/load_later/homopolymer_analysis/all_methods_info_mtFilteredV5_filterVariantDensityIsoSeq_noNearSJ.RData")
# save(dat_full, file="~/load_later/homopolymer_analysis/all_methods_info_mtFilteredV5_filterVariantDensityIsoSeq_noNearSJ_isoSeqCover20.RData")


