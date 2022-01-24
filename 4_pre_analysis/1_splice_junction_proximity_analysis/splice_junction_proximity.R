# this script is to create and save the data needed to make the figures  to analyze the effect of 
# splice junction proximity on variant calling performance


### load packages
library(lrRNAseqBenchmark)


### for snps
master_table1 <- "/home/vbarbo/project_2021/paper_analysis/extra_files/master_tables/mt_jurkat_allMethods_filtered_v7.RData"
master_table2 <- "/home/vbarbo/project_2021/paper_analysis/extra_files/master_tables/mt_wtc11_allMethods_filtered_v7.RData"
min_isoseq_coverage <- 20
method_names <- c("dv_s_fc", "c3_mix", "gatk_s")
variant_type <- "snp"
output_method_names <- c("SNCR+FC+DeepVariant", "Clair3-mix", "SNCR+GATK")
truth_names <- c("merged", "allen")
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
variant_type <- "indel"

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
     file="/home/vbarbo/project_2021/paper_analysis/extra_files/splice_junction_proximity_analysis.Rdata")

