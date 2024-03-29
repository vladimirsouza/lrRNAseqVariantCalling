# this script is to create and save the data needed to make the figures  to analyze the effect of 
# splice junction proximity on variant calling performance


### load packages
library(lrRNAseqBenchmark)


### for snps
rds_master_table1 <- "/home/vbarbo/project_2021/paper_analysis/extra_files/master_tables/mt_jurkat_allMethods_filtered_v7.rds"
rds_master_table2 <- "/home/vbarbo/project_2021/paper_analysis/extra_files/master_tables/mt_wtc11_allMethods_filtered_v7.rds"
min_isoseq_coverage <- 20
method_names <- c("dv_s_fc", "c3_mix", "gatk_s")
variant_type <- "snp"
output_method_names <- c("SNCR+FC+DeepVariant", "Clair3-mix", "SNCR+GATK")
truth_names <- c("merged", "allen")
experiment_names <- c("Jurkat", "WTC-11")
method_dataset_name <- c("isoSeq", "isoSeq")


master_table1 <- readRDS(rds_master_table1)
master_table2 <- readRDS(rds_master_table2)

sj_proximity_snps <- splice_junction_analysis_table(master_table1,
                                                    master_table2,
                                                    experiment_names=experiment_names,
                                                    truth_names=truth_names,
                                                    method_dataset_name=method_dataset_name,
                                                    method_names=method_names,
                                                    output_method_names=output_method_names,
                                                    variant_type=variant_type,
                                                    min_isoseq_coverage=min_isoseq_coverage)

### for indels
variant_type <- "indel"

sj_proximity_indels <- splice_junction_analysis_table(master_table1,
                                                      master_table2,
                                                      experiment_names=experiment_names,
                                                      truth_names=truth_names,
                                                      method_dataset_name=method_dataset_name,
                                                      method_names=method_names,
                                                      output_method_names=output_method_names,
                                                      variant_type=variant_type,
                                                      min_isoseq_coverage=min_isoseq_coverage)

### save the objects to a file
saveRDS(sj_proximity_snps, "/home/vbarbo/project_2021/paper_analysis/extra_files/snp_splice_junction_proximity_analysis.rds")
saveRDS(sj_proximity_indels, "/home/vbarbo/project_2021/paper_analysis/extra_files/indel_splice_junction_proximity_analysis.rds")

