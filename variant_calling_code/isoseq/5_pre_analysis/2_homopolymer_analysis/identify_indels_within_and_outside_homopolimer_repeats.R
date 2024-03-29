# this script is to create the table that stores the data needed to draw  the plot for the homopolymer analysis

# we kept only sites that:
#   * commun filter:
#     * short-read variant density: x <= 3                          (ignore positions)
#     * ignore intronic regions (isoSeq_coverage > 0)               (ignore positions)
#     * DeepVariant QUAL: x >= 15                                   (remove method calling)
#     * short-read coverare : quantil 5% <= x <= 95%                (ignore positions)
#   * extra filtering:
#     * iso-seq variant density: x <= 3
#     * are far from splice junction, i.e. further than 20 bp;
#     * iso-seq coverage is higher than 20 reads (by function `method_homopolymer_indels`);
#     * are classified as heterozygous alternative (by function `method_homopolymer_indels`);
#     * indels (by function `method_homopolymer_indels`);


### load packages
library(lrRNAseqBenchmark)
library(dplyr)
library(Biostrings)
library(snakecase)


### IMPORTANT: there are other variables defined through the script
# METHOD_NAMES <- c("dv_s_fc", "c3_mix", "gatk_s")
METHOD_VCF_FILES <- c(
  "/home/vbarbo/project_2021/paper_analysis/wtc11/variant_calling_from_isoseq/deepvariant/dv_sncr_fc/deepvariant_calls_pass.vcf.gz",
  "/home/vbarbo/project_2021/paper_analysis/wtc11/variant_calling_from_isoseq/clair3/mix/pileup_pass_mix_norep.recode.vcf.gz",
  "/home/vbarbo/project_2021/paper_analysis/wtc11/variant_calling_from_isoseq/gatk/isoSeq_wtc11.recal_pass.vcf.gz"
)
METHOD_DATASET_NAME <- "isoSeq"
TRUTH_NAME <- "allen"
TRUTH_VCF_FILE <- "/home/vbarbo/project_2021/paper_analysis/wtc11/ground_truth/3546dc62_AH77TTBBXX_DS-229105_GCCAAT_recalibrated_subsetChromosomes_pass.vcf.gz"
MASTER_TABLE <- "/home/vbarbo/project_2021/paper_analysis/extra_files/master_tables/mt_wtc11_allMethods_filtered_v7.rds"
REF_FASTA <- "/home/vbarbo/project_2021/paper_analysis/reference/genome/GRCh38.p13_all_chr.fasta"



### load the filtered master table and add additional filtering:
# * filter out sites in high-variant-density regions of Iso-Seq read alignemnts;
# * filter out sites near splice junctions.

dat1 <- readRDS(MASTER_TABLE)
dim(dat1)

dat1 <- filter(dat1, variantDensity_dv_dvS_dvSFc_c3_c3S_c3SFc_c3Mix_gatkS_ncS <= 3)
dim(dat1)

dat1 <- filter(dat1, is_near_ss==0)
dim(dat1)



### load reference genome sequences
ref_fasta_seqs <- readDNAStringSet(REF_FASTA)
names(ref_fasta_seqs) <- sub("(^chr[0-9]+|X|Y).*", "\\1", names(ref_fasta_seqs))

# # get all homopolymers of the reference genome
# homopolymers <- homopolymerFinder(ref_fasta_seqs)
### load file with all homopolymers of the reference genome created in script /home/vbarbo/project_2021/projects/lrRNAseqVariantCalling/4_create_master_table/jurkat/mt_jurkat_allMethods_v7.Rmd
homopolymers <- readRDS("/home/vbarbo/project_2021/paper_analysis/extra_files/homopolymers_GRCh38.p13_all_chr.rds")



### get indels in homopolymers (homopolymer length equal to 1 means non homopolymer)
### also, filter sites by minimum iso-seq read coverage, keep only heterozygous indels

# dv_s_fc
dat_hom_dvSFc <- method_homopolymer_indels(input_table=dat1,
                                           first_method_name=TRUTH_NAME,
                                           second_method_name="dv_s_fc",
                                           vcf_first=TRUTH_VCF_FILE,
                                           vcf_second=METHOD_VCF_FILES[1],
                                           method_dataset_name=METHOD_DATASET_NAME,
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
                                           method_dataset_name=METHOD_DATASET_NAME,
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
                                           method_dataset_name=METHOD_DATASET_NAME,
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
  input_hom_table = dat_hom_dvSFc,
  variant_type = "deletion",
  method_name = "dv_s_fc",
  truth_name = "allen",
  hom_length_intervals = c(1,2,3,5,11,16,21),
  interval_names = c("non-hp",2,"3-4","5-10", "11-15", "16-20", ">=21"),
  to_calculate = "pre_rec_f1",
  output_method_name = "SNCR+FC+DeepVariant"
)
dat_full <- mapply(rbind, dat_full, k, SIMPLIFY=FALSE)

# insertion
k <- make_homopolymer_table_to_plot(
  input_hom_table = dat_hom_dvSFc,
  variant_type = "insertion",
  method_name = "dv_s_fc",
  truth_name = "allen",
  hom_length_intervals = c(1,2,3,5,11,16,21),
  interval_names = c("non-hp",2,"3-4","5-10", "11-15", "16-20", ">=21"),
  to_calculate = "pre_rec_f1",
  output_method_name = "SNCR+FC+DeepVariant"
)
dat_full <- mapply(rbind, dat_full, k, SIMPLIFY=FALSE)

# # # # # # # # # #
# # clair3 mix  # #
# # # # # # # # # #

k <- make_homopolymer_table_to_plot(
  input_hom_table = dat_hom_c3Mix,
  variant_type = "deletion",
  method_name = "c3_mix",
  truth_name = "allen",
  hom_length_intervals = c(1,2,3,5,11,16,21),
  interval_names = c("non-hp",2,"3-4","5-10", "11-15", "16-20", ">=21"),
  to_calculate = "pre_rec_f1",
  output_method_name = "Clair3 mix"
)
dat_full <- mapply(rbind, dat_full, k, SIMPLIFY=FALSE)

k <- make_homopolymer_table_to_plot(
  input_hom_table = dat_hom_c3Mix,
  variant_type = "insertion",
  method_name = "c3_mix",
  truth_name = "allen",
  hom_length_intervals = c(1,2,3,5,11,16,21),
  interval_names = c("non-hp",2,"3-4","5-10", "11-15", "16-20", ">=21"),
  to_calculate = "pre_rec_f1",
  output_method_name = "Clair3 mix"
)
dat_full <- mapply(rbind, dat_full, k, SIMPLIFY=FALSE)

# # # # # # # # #
# # sncr+gatk # #
# # # # # # # # #

k <- make_homopolymer_table_to_plot(
  input_hom_table = dat_hom_gatkS,
  variant_type = "deletion",
  method_name = "gatk_s",
  truth_name = "allen",
  hom_length_intervals = c(1,2,3,5,11,16,21),
  interval_names = c("non-hp",2,"3-4","5-10", "11-15", "16-20", ">=21"),
  to_calculate = "pre_rec_f1",
  output_method_name = "SNCR+GATK"
)
dat_full <- mapply(rbind, dat_full, k, SIMPLIFY=FALSE)

k <- make_homopolymer_table_to_plot(
  input_hom_table = dat_hom_gatkS,
  variant_type = "insertion",
  method_name = "gatk_s",
  truth_name = "allen",
  hom_length_intervals = c(1,2,3,5,11,16,21),
  interval_names = c("non-hp",2,"3-4","5-10", "11-15", "16-20", ">=21"),
  to_calculate = "pre_rec_f1",
  output_method_name = "SNCR+GATK"
)
dat_full <- mapply(rbind, dat_full, k, SIMPLIFY=FALSE)



### some manipulations
k <- recode(levels(dat_full$class_counts$Measures),
            "precision"="Precision", "sensitivity"="Recall", "f1Score"="F1-score")
levels(dat_full$class_counts$Measures) <- k

dat_full$dat_text$y <- max(dat_full$dat_text$y)

k <- c("SNCR+FC+DeepVariant", "Clair3 mix", "SNCR+GATK")
dat_full$class_counts$method <- factor(dat_full$class_counts$method,
                                       levels=k, ordered=TRUE)
dat_full$dat_text$method <- factor(dat_full$dat_text$method,
                                   levels=k, ordered=TRUE)


dat_homopolymer <- dat_full

### save object `dat_homopolymer` to a file.
### this is the data used to draw the chart for the homopolymer analysis
saveRDS(dat_homopolymer, "/home/vbarbo/project_2021/paper_analysis/extra_files/dat_homopolymer_analysis.rds")
