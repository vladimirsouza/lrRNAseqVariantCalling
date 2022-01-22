# this script is to create and save the data needed to make the figure for the analysis of
# the effect of splice junction proximity on varaint calling performance


### load packages
library(variantCallingFromIsoSeq)


### for snps
master_table1 <- "/home/vbarbo/project_2021/projects/lrRNA-seq_variant_calling/3_creating_master_table/jurkat/mt_jurkat_allMethods_filtered_v6.RData"
master_table2 <- "/home/vbarbo/project_2021/projects/lrRNA-seq_variant_calling/3_creating_master_table/wtc11/mt_wtc11_allMethods_filtered_v6.RData"
min_isoseq_coverage <- 20
method_names <- c("dv_s_fc", "c3_mix", "gatk_s")
variant_type <- "snp"
output_method_names <- c("SNCR+FC+DeepVariant", "Clair3 mix", "SNCR+GATK")
truth_names <- c("jurkat_dna_merged", "allen")
experiment_names <- c("Jurkat", "WTC-11")

sj_proximity_snps <- splice_junction_analysis_table(master_table1,
                                                    master_table2,
                                                    experiment_names=experiment_names,
                                                    truth_names=truth_names,
                                                    method_names=method_names,
                                                    output_method_names=output_method_names,
                                                    variant_type=variant_type,
                                                    min_isoseq_coverage=min_isoseq_coverage)

### for indels
master_table1 <- "~/project_2021/scripts/paper_1/master_tables/jurkat/mt_jurkat_allMethods_filtered_v5.RData"
master_table2 <- "~/project_2021/scripts/paper_1/master_tables/wtc11/mt_wtc11_allMethods_filtered_v5.RData"
min_isoseq_coverage <- 20
method_names <- c("dv_s_fc", "c3_mix", "gatk_s")
variant_type <- "indel"
output_method_names <- c("SNCR+FC+DeepVariant", "Clair3 mix", "SNCR+GATK")
truth_names <- c("jurkat_dna_merged", "allen")
experiment_names <- c("Jurkat", "WTC-11")

sj_proximity_indels <- splice_junction_analysis_table(master_table1,
                                                     master_table2,
                                                     experiment_names=experiment_names,
                                                     truth_names=truth_names,
                                                     method_names=method_names,
                                                     output_method_names=output_method_names,
                                                     variant_type=variant_type,
                                                     min_isoseq_coverage=min_isoseq_coverage)

### save the objects to a file
save(sj_proximity_snps, sj_proximity_indels,
     file="~/load_later/near_splice_junctions/data_needed_to_make_the_chart.Rdata")

