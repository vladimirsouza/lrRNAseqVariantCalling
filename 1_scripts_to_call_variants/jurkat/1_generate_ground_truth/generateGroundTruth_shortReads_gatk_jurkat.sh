#!/bin/bash


### this script is to generate our ground-truth vcf file from short-read data.
### we download the illumina data from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5941560/

### in this same paper they also
### generate the vcf file, but it's tetraploid. but they provide the code that
### generates the vcf, so i can get advantage of it to create my own diploid
### ground-truth vcf file.

### good documentation in:
### A guide to GATK4 best practice pipeline performance and optimization on the IBM OpenPOWER system
### https://www.ibm.com/downloads/cas/ZJQD0QAL
### also -- paper from which I downloaded the dna short-read data to create the ground truth:
### https://bitbucket.org/sulab/jurkat_variant_calling/src/master/




#################################################################################
###### get the short-read fatsq files ###########################################
###### there are two datasets: jurkat_wgs_pe_100bp and jurkat_wgs_pe_150bp ######
#################################################################################


### downloading the datasets (there are two differente illumina data files for jurkat cells)
cd /home/vbarbo/project_2021/datasets/gloria_data/dna_short_read_jurkat_downloaded/jurkat_wgs_pe_100bp
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-7/SRR5349450/SRR5349450.1
mv SRR5349450.1 jurkat_wgs_pe_100bp

cd /home/vbarbo/project_2021/datasets/gloria_data/dna_short_read_jurkat_downloaded/jurkat_wgs_pe_150bp
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-7/SRR5349449/SRR5349449.1
mv SRR5349449.1 jurkat_wgs_pe_150bp



### Convert SRA data into fastq format
cd /home/vbarbo/project_2021/datasets/gloria_data/dna_short_read_jurkat_downloaded/jurkat_wgs_pe_100bp
fastq-dump -I --split-files jurkat_wgs_pe_100bp
cd /home/vbarbo/project_2021/datasets/gloria_data/dna_short_read_jurkat_downloaded/jurkat_wgs_pe_150bp
fastq-dump -I --split-files jurkat_wgs_pe_150bp


### I shouldn't have used the fastq-dump argument -I. Now I have to remove the last number in fastq sequence name (before first space character)
cd /home/vbarbo/project_2021/datasets/gloria_data/dna_short_read_jurkat_downloaded/jurkat_wgs_pe_100bp
sed -i '/^[@+]jurkat_wgs_pe_100bp/ s/\.1 / /' jurkat_wgs_pe_100bp_1.fastq
sed -i '/^[@+]jurkat_wgs_pe_100bp/ s/\.2 / /' jurkat_wgs_pe_100bp_2.fastq
cd /home/vbarbo/project_2021/datasets/gloria_data/dna_short_read_jurkat_downloaded/jurkat_wgs_pe_150bp
sed -i '/^[@+]jurkat_wgs_pe_150bp/ s/\.1 / /' jurkat_wgs_pe_150bp_1.fastq
sed -i '/^[@+]jurkat_wgs_pe_150bp/ s/\.2 / /' jurkat_wgs_pe_150bp_2.fastq






##########################################
###### jurkat_wgs_pe_100bp ###############
###### align reads to the reference ######
##########################################


### inputs
READ1=/home/vbarbo/project_2021/datasets/gloria_data/dna_short_read_jurkat_downloaded/jurkat_wgs_pe_100bp/jurkat_wgs_pe_100bp_1.fastq
READ2=/home/vbarbo/project_2021/datasets/gloria_data/dna_short_read_jurkat_downloaded/jurkat_wgs_pe_100bp/jurkat_wgs_pe_100bp_2.fastq
REF=/home/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta
OUTPUT_DIR=/home/vbarbo/project_2021/datasets/gloria_data/analysis/my_ground_truth_jurkat_wgs_pe_100bp
SAMPLE=jurkat100bp
THREADS=30

#DBSNP=/home/vbarbo/project_2021/datasets/reference/dbSNP_GRCh38/00-common_all.vcf.gz
#DBSNP=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/Homo_sapiens_assembly38.dbsnp138.vcf.gz
DBSNP=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/Homo_sapiens_assembly38.known_indels.vcf.gz
OMNI=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/1000G_omni2.5.hg38.vcf.gz
G1000=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz
HAPMAP=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/hapmap_3.3.hg38.vcf.gz
MILLS=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz



### index the ref
bwa index $REF


### align reads to the reference
filename=$SAMPLE
bwa mem -t $THREADS -Ma \
  -R @RG\\tID:${filename}\\tSM:${filename}\\tPL:ILM\\tLB:${filename} \
  $REF \
  $READ1 \
  $READ2 \
  > ${OUTPUT_DIR}/${SAMPLE}_aln.sam
#  | samtools sort - -@ $THREADS -n -m 4G -o ${OUTPUT_DIR}/${SAMPLE}_sorted.bam


### create read group
######### <<<<<<<<<<<<<<<-----------------------------=========== can i use multiple cores here?
#java -XX:ParallelGCThreads=$THREADS -jar /home/vbarbo/programs/picard/build/libs/picard.jar AddOrReplaceReadGroups \
java -jar /home/vbarbo/programs/picard/build/libs/picard.jar AddOrReplaceReadGroups \
  -I ${OUTPUT_DIR}/${SAMPLE}_aln.sam \
  -O ${OUTPUT_DIR}/${SAMPLE}_sorted.bam \
  -SORT_ORDER queryname \
  -VALIDATION_STRINGENCY LENIENT \
  -MAX_RECORDS_IN_RAM 5000000 \
  -ID 1 \
  -LB CTTGTA \
  -PL ILLUMINA \
  -PU 1234 \
  -SM jurkat
#  -SORT_ORDER coordinate


### Sequence-marking duplicates
samtools fixmate -m -@ $THREADS \
  ${OUTPUT_DIR}/${SAMPLE}_sorted.bam \
  ${OUTPUT_DIR}/${SAMPLE}_fixmate.bam
samtools sort -@ $THREADS -m 6G \
  -o ${OUTPUT_DIR}/${SAMPLE}_sorted.bam \
  ${OUTPUT_DIR}/${SAMPLE}_fixmate.bam
samtools markdup -s -@ $THREADS \
  ${OUTPUT_DIR}/${SAMPLE}_sorted.bam \
  ${OUTPUT_DIR}/${SAMPLE}_sorted_dedup.bam
samtools index -@ $THREADS \
  ${OUTPUT_DIR}/${SAMPLE}_sorted_dedup.bam


# ## picard mark duplicates (is this marking all of them with `PG:Z:MarkDuplicates` ?)
# java -jar /home/vbarbo/programs/picard/build/libs/picard.jar MarkDuplicates \
#   -I ${OUTPUT_DIR}/${SAMPLE}_sorted.bam \
#   -O ${OUTPUT_DIR}/${SAMPLE}_dedup.bam \
#   -M ${OUTPUT_DIR}/${SAMPLE}_metrics.txt \
#   -ASSUME_SORTED true \
#   -VALIDATION_STRINGENCY LENIENT


# ## index the bam file
# java -jar /home/vbarbo/programs/picard/build/libs/picard.jar BuildBamIndex \
#   -I ${OUTPUT_DIR}/${SAMPLE}_dedup.bam






##########################################
###### jurkat_wgs_pe_150bp ###############
###### align reads to the reference ######
##########################################


### inputs
READ1=/home/vbarbo/project_2021/datasets/gloria_data/dna_short_read_jurkat_downloaded/jurkat_wgs_pe_150bp/jurkat_wgs_pe_150bp_1.fastq
READ2=/home/vbarbo/project_2021/datasets/gloria_data/dna_short_read_jurkat_downloaded/jurkat_wgs_pe_150bp/jurkat_wgs_pe_150bp_2.fastq
REF=/home/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta
OUTPUT_DIR=/home/vbarbo/project_2021/datasets/gloria_data/analysis/my_ground_truth_jurkat_wgs_pe_150bp
SAMPLE=jurkat150bp
THREADS=30

#DBSNP=/home/vbarbo/project_2021/datasets/reference/dbSNP_GRCh38/00-common_all.vcf.gz
#DBSNP=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/Homo_sapiens_assembly38.dbsnp138.vcf.gz
DBSNP=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/Homo_sapiens_assembly38.known_indels.vcf.gz
OMNI=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/1000G_omni2.5.hg38.vcf.gz
G1000=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz
HAPMAP=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/hapmap_3.3.hg38.vcf.gz
MILLS=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz



### index the ref
bwa index $REF


### align reads to the reference
filename=$SAMPLE
bwa mem -t $THREADS -Ma \
  -R @RG\\tID:${filename}\\tSM:${filename}\\tPL:ILM\\tLB:${filename} \
  $REF \
  $READ1 \
  $READ2 \
  > ${OUTPUT_DIR}/${SAMPLE}_aln.sam
#  | samtools sort - -@ $THREADS -n -m 4G -o ${OUTPUT_DIR}/${SAMPLE}_sorted.bam


### create read group
######### <<<<<<<<<<<<<<<-----------------------------=========== can i use multiple cores here?
#java -XX:ParallelGCThreads=$THREADS -jar /home/vbarbo/programs/picard/build/libs/picard.jar AddOrReplaceReadGroups \
java -jar /home/vbarbo/programs/picard/build/libs/picard.jar AddOrReplaceReadGroups \
  -I ${OUTPUT_DIR}/${SAMPLE}_aln.sam \
  -O ${OUTPUT_DIR}/${SAMPLE}_sorted.bam \
  -SORT_ORDER queryname \
  -VALIDATION_STRINGENCY LENIENT \
  -MAX_RECORDS_IN_RAM 5000000 \
  -ID 1 \
  -LB CTTGTA \
  -PL ILLUMINA \
  -PU 1234 \
  -SM jurkat
#  -SORT_ORDER coordinate


### Sequence-marking duplicates
samtools fixmate -m -@ $THREADS \
  ${OUTPUT_DIR}/${SAMPLE}_sorted.bam \
  ${OUTPUT_DIR}/${SAMPLE}_fixmate.bam
samtools sort -@ $THREADS -m 6G \
  -o ${OUTPUT_DIR}/${SAMPLE}_sorted.bam \
  ${OUTPUT_DIR}/${SAMPLE}_fixmate.bam
samtools markdup -s -@ $THREADS \
  ${OUTPUT_DIR}/${SAMPLE}_sorted.bam \
  ${OUTPUT_DIR}/${SAMPLE}_sorted_dedup.bam
samtools index -@ $THREADS \
  ${OUTPUT_DIR}/${SAMPLE}_sorted_dedup.bam






#########################################################################################
###### merge dna_100 and dna_150 bam files to call variants with high(er) coverage ######
#########################################################################################


### inputs
BAM1=/home/vbarbo/project_2021/datasets/gloria_data/analysis/my_ground_truth_jurkat_wgs_pe_100bp/jurkat100bp_sorted_dedup.bam
BAM2=/home/vbarbo/project_2021/datasets/gloria_data/analysis/my_ground_truth_jurkat_wgs_pe_150bp/jurkat150bp_sorted_dedup.bam
REF=/home/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta
OUTPUT_DIR=/home/vbarbo/project_2021/datasets/gloria_data/analysis/truth_merged_bams
SAMPLE=jurkat150bp # i forgot to change the sample name to merged.
                   # but keep in mind that ~/project_2021/datasets/gloria_data/analysis/truth_merged_bams/primary/jurkat150bp.recal_pass.vcf
                   # contains not only varaints from jurkat150bp, but also from jurkat100bp
THREADS=30

DBSNP=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/Homo_sapiens_assembly38.known_indels.vcf.gz
OMNI=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/1000G_omni2.5.hg38.vcf.gz
G1000=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz
HAPMAP=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/hapmap_3.3.hg38.vcf.gz
MILLS=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz




### merge bams
samtools merge \
  -@ $THREADS \
  ${OUTPUT_DIR}/${SAMPLE}.bam \
  $BAM1 \
  $BAM2


### remove non-primary aligments (i.e., supplementary and secundary)
samtools view -S -b -F 2308 -@ $THREADS \
  -o ${OUTPUT_DIR}/primary/${SAMPLE}_primary.bam \
  ${OUTPUT_DIR}/${SAMPLE}.bam


### index bam
samtools index -@ $THREADS \
  ${OUTPUT_DIR}/primary/${SAMPLE}_primary.bam







###############################################
###### call variants from the merged bam ######
###############################################

### more inputs
OUTPUT_DIR=/home/vbarbo/project_2021/datasets/gloria_data/analysis/truth_merged_bams/primary


# https://www.biostars.org/p/363063/
# Indel Realignment is no longer neccassery and recommended by GATK, if you use the HaplotypeCaller.

# https://github.com/broadinstitute/gatk-docs/blob/master/blog-2012-to-2019/2016-06-21-Changing_workflows_around_calling_SNPs_and_indels.md?id=7847
# "As announced in the GATK v3.6 highlights, variant calling workflows that use HaplotypeCaller or MuTect2 now omit
#  indel realignment."
# "Realigning reads using IndelRealigner or assembling reads using HaplotypeCaller allows us to call the insertion."
# "HaplotypeCaller reassembles and calls this insertion with or without indel realignment."
# "If indel realignment is redundant to HaplotypeCaller’s functions, why were we recommending both until now? Well,
#  although HaplotypeCaller throws out realignment’s local placement of reads, indel realignment still influences
#  outcomes. I’m told these influences should be subtle for high quality data but make a difference for lower
#  quality data."

# https://www.biostars.org/p/305123/
# "Basically, the realignment step has been integrated into the GATK variant callers and they can output their own
#  realigned BAMs. However, regardless of how you perform realignment, it is still preferred."

# gatk --java-options "-Xmx4g" BaseRecalibrator -h

# https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator
# "First pass of the base quality score recalibration."
# "Generates a recalibration table based on various covariates."
# "It does a by-locus traversal operating only at sites that are in the known sites VCF. ExAc, gnomAD, or dbSNP
#  resources can be used as known sites of variation."
# "We assume that all reference mismatches we see are therefore errors and indicative of poor base quality. Since
#  there is a large amount of data one can then calculate an empirical probability of error given the particular 
#  covariates seen at this site, where p(error) = num mismatches / num observations".
# "The output file is a table"

# https://www.intel.com/content/dam/www/public/us/en/documents/white-papers/deploying-gatk-best-practices-paper.pdf
# [-L exome_targets.intervals] #if running on exome
# [-ip 50] #if running on exome (interval padding)





### to use multiple cores, create intervals
interval_root_dir=/home/vbarbo/project_2021/datasets/gloria_data/analysis/my_ground_truth_jurkat_wgs_pe_100bp
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" ScatterIntervalsByNs \
  -R $REF \
  -O $interval_root_dir/jurkat100bp.interval_list
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" SplitIntervals \
  -R $REF \
  -L $interval_root_dir/jurkat100bp.interval_list \
  --scatter-count $THREADS \
  -O $interval_root_dir/jurkat100bp_scattered.interval_list

SCATTERED_INTERVAL_LIST=/home/vbarbo/project_2021/datasets/gloria_data/analysis/my_ground_truth_jurkat_wgs_pe_100bp/jurkat100bp_scattered.interval_list
THREADS=30
loop_num=`expr $THREADS - 1`
# as long as we are using the same /home/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta as the reference
# and 30 cores, it's okay to use this $SCATTERED_INTERVAL_LIST for other analysis



### calculate base recalibration
mkdir ${OUTPUT_DIR}/base_recalibration_tables
for i in `seq -f '%04g' 0 $loop_num`
do
  gatk --java-options "-Xmx4G" BaseRecalibrator \
    -R $REF \
    -I ${OUTPUT_DIR}/${SAMPLE}.bam \
    -O ${OUTPUT_DIR}/base_recalibration_tables/${SAMPLE}_recal_data_$i.table \
    -L $SCATTERED_INTERVAL_LIST/$i-scattered.interval_list \
    --known-sites $DBSNP \
    --known-sites $MILLS \
    --known-sites $G1000 &
######## <<<<<<-----------==== first, try without `--maximum-cycle-value 100000`
#    --maximum-cycle-value 100000 \
#    -L ${OUTPUT_DIR}/${SAMPLE}_scattered.interval_list/$i-scattered.interval_list \
done
wait


### apply base recalibration
mkdir ${OUTPUT_DIR}/base_recalibration_bams
for i in `seq -f '%04g' 0 $loop_num`
do
  gatk --java-options "-Xmx4G" ApplyBQSR \
  -R $REF \
  -I ${OUTPUT_DIR}/${SAMPLE}.bam \
  -bqsr ${OUTPUT_DIR}/base_recalibration_tables/${SAMPLE}_recal_data_$i.table \
  -O ${OUTPUT_DIR}/base_recalibration_bams/${SAMPLE}_recal_$i.bam \
  -L $SCATTERED_INTERVAL_LIST/$i-scattered.interval_list \
  --static-quantized-quals 10 \
  --static-quantized-quals 20 \
  --static-quantized-quals 30 &
#  -L ${OUTPUT_DIR}/${SAMPLE}_scattered.interval_list/$i-scattered.interval_list \
done
wait



### Variant calling using GATK HaplotypeCaller with GVCF mode
mkdir ${OUTPUT_DIR}/intermediate_gvcfs
for i in `seq -f '%04g' 0 $loop_num`
do
  gatk --java-options "-Xmx4G" HaplotypeCaller \
    -R $REF \
    -I ${OUTPUT_DIR}/base_recalibration_bams/${SAMPLE}_recal_$i.bam \
    -O ${OUTPUT_DIR}/intermediate_gvcfs/${SAMPLE}_recal_$i.g.vcf \
    -L $SCATTERED_INTERVAL_LIST/$i-scattered.interval_list \
    --native-pair-hmm-threads 1 \
    -ERC GVCF &
#    -L ${OUTPUT_DIR}/${SAMPLE}_scattered.interval_list/$i-scattered.interval_list \
#    -pairHMM VSX_LOGLESS_CACHING
#    -stand-call-conf 10
done
wait


### Consolidate and genotype genomic variant call formats (GVCFs)
# genotype gvcf files
mkdir ${OUTPUT_DIR}/vcf_parts
for i in `seq -f '%04g' 0 $loop_num`
do
  gatk --java-options "-Xmx4G" GenotypeGVCFs \
  -R $REF \
  -V ${OUTPUT_DIR}/intermediate_gvcfs/${SAMPLE}_recal_$i.g.vcf \
  -L $SCATTERED_INTERVAL_LIST/$i-scattered.interval_list \
  -O ${OUTPUT_DIR}/vcf_parts/${SAMPLE}_variants_$i.vcf &
#  -L ${OUTPUT_DIR}/${SAMPLE}_scattered.interval_list/$i-scattered.interval_list \
done
wait


# merge scattered phenotype vcf files
vcfFilesArg=(${OUTPUT_DIR}/vcf_parts/${SAMPLE}_variants_*.vcf)
vcfFilesArg=${vcfFilesArg[@]/#/-I }
gatk --java-options "-Xmx4G" GatherVcfs \
  -R $REF \
  $vcfFilesArg \
  -O ${OUTPUT_DIR}/${SAMPLE}.vcf


### VARIANT QUALITY SCORE RECALIBRATION - SNPs
mkdir ${OUTPUT_DIR}/variant_recalibration
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" VariantRecalibrator \
  -V ${OUTPUT_DIR}/${SAMPLE}.vcf \
  -O ${OUTPUT_DIR}/variant_recalibration/${SAMPLE}_recalibrate_SNP.recal \
  -mode SNP \
  --tranches-file ${OUTPUT_DIR}/variant_recalibration/${SAMPLE}_recalibrate_SNP.tranches \
  -tranche 100.0 \
  -tranche 99.9 \
  -tranche 99.0 \
  -tranche 90.0 \
  -an QD \
  -an FS \
  -an MQRankSum \
  -an ReadPosRankSum \
  -an SOR \
  -an MQ \
  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
  -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI \
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1000 \
  -resource:dbsnp,known=true,training=false,truth=false,prior=7.0 $DBSNP
#  --max-gaussians 6


### Apply recalibration to SNPs
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" ApplyVQSR \
  -V ${OUTPUT_DIR}/${SAMPLE}.vcf \
  -O ${OUTPUT_DIR}/variant_recalibration/${SAMPLE}_recalibrated_snps_raw_indels.vcf \
  --recal-file ${OUTPUT_DIR}/variant_recalibration/${SAMPLE}_recalibrate_SNP.recal \
  --tranches-file ${OUTPUT_DIR}/variant_recalibration/${SAMPLE}_recalibrate_SNP.tranches \
  -truth-sensitivity-filter-level 99.5 \
  --create-output-variant-index true \
  -mode SNP


### Run Variant Recalibrator – Indels
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" VariantRecalibrator \
  -V ${OUTPUT_DIR}/variant_recalibration/${SAMPLE}_recalibrated_snps_raw_indels.vcf \
  -O ${OUTPUT_DIR}/variant_recalibration/${SAMPLE}_recalibrate_INDEL.recal \
  -mode INDEL \
  --tranches-file ${OUTPUT_DIR}/variant_recalibration/${SAMPLE}_recalibrate_INDEL.tranches \
  -tranche 100.0 \
  -tranche 99.9 \
  -tranche 99.0 \
  -tranche 90.0 \
  -an QD \
  -an FS \
  -an MQRankSum \
  -an ReadPosRankSum \
  -an SOR \
  -resource:mills,known=false,training=true,truth=true,prior=12.0 $MILLS \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP
#  --max-gaussians 4


### Apply recalibration to Indels
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" ApplyVQSR \
  -V ${OUTPUT_DIR}/variant_recalibration/${SAMPLE}_recalibrated_snps_raw_indels.vcf \
  -O ${OUTPUT_DIR}/${SAMPLE}.recal.vcf \
  --recal-file ${OUTPUT_DIR}/variant_recalibration/${SAMPLE}_recalibrate_INDEL.recal \
  --tranches-file ${OUTPUT_DIR}/variant_recalibration/${SAMPLE}_recalibrate_INDEL.tranches \
  -truth-sensitivity-filter-level 99.0 \
  --create-output-variant-index true \
  -mode INDEL

