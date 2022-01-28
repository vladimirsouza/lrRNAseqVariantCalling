#!/bin/bash


### this script is to generate the ground-truth vcf file (Jurkat dataset) from short-read data.
### we download illumina data from the paper in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5941560/

### in this same paper they also
### generate the vcf file, but it's tetraploid. but they provide the code that
### generates the vcf, so we can get advantage of it to create our own diploid
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
cd /home/vbarbo/project_2021/paper_analysis/jurkat/data/dna_short_reads/jurkat_wgs_pe_100bp
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-7/SRR5349450/SRR5349450.1
mv SRR5349450.1 jurkat_wgs_pe_100bp

cd /home/vbarbo/project_2021/paper_analysis/jurkat/data/dna_short_reads/jurkat_wgs_pe_150bp
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-7/SRR5349449/SRR5349449.1
mv SRR5349449.1 jurkat_wgs_pe_150bp



### Convert SRA data into fastq format
cd /home/vbarbo/project_2021/paper_analysis/jurkat/data/dna_short_reads/jurkat_wgs_pe_100bp
fastq-dump --split-files jurkat_wgs_pe_100bp
cd /home/vbarbo/project_2021/paper_analysis/jurkat/data/dna_short_reads/jurkat_wgs_pe_150bp
fastq-dump --split-files jurkat_wgs_pe_150bp

### compress files
bgzip /home/vbarbo/project_2021/paper_analysis/jurkat/data/dna_short_reads/jurkat_wgs_pe_100bp/jurkat_wgs_pe_100bp_1.fastq
bgzip /home/vbarbo/project_2021/paper_analysis/jurkat/data/dna_short_reads/jurkat_wgs_pe_100bp/jurkat_wgs_pe_100bp_2.fastq

bgzip /home/vbarbo/project_2021/paper_analysis/jurkat/data/dna_short_reads/jurkat_wgs_pe_150bp/jurkat_wgs_pe_150bp_1.fastq
bgzip /home/vbarbo/project_2021/paper_analysis/jurkat/data/dna_short_reads/jurkat_wgs_pe_150bp/jurkat_wgs_pe_150bp_2.fastq




##########################################
###### jurkat_wgs_pe_100bp ###############
###### align reads to the reference ######
##########################################


### inputs
READ1=/home/vbarbo/project_2021/paper_analysis/jurkat/data/dna_short_reads/jurkat_wgs_pe_100bp/jurkat_wgs_pe_100bp_1.fastq.gz
READ2=/home/vbarbo/project_2021/paper_analysis/jurkat/data/dna_short_reads/jurkat_wgs_pe_100bp/jurkat_wgs_pe_100bp_2.fastq.gz
REF=/home/vbarbo/project_2021/paper_analysis/reference/genome/GRCh38.p13_all_chr.fasta
OUTPUT_DIR=/home/vbarbo/project_2021/paper_analysis/jurkat/ground_truth
SAMPLE1=jurkat100bp
THREADS=30




### index the ref
bwa index $REF

samtools faidx $REF

gatk CreateSequenceDictionary -R $REF


### align reads to the reference
bwa mem -t $THREADS -Ma \
  -R @RG\\tID:${SAMPLE1}\\tSM:${SAMPLE1}\\tPL:ILM\\tLB:${SAMPLE1} \
  $REF \
  $READ1 \
  $READ2 \
  > ${OUTPUT_DIR}/${SAMPLE1}_aln.sam


### create read group
java -XX:ParallelGCThreads=$THREADS -jar /home/vbarbo/programs/picard/build/libs/picard.jar AddOrReplaceReadGroups \
  -I ${OUTPUT_DIR}/${SAMPLE1}_aln.sam \
  -O ${OUTPUT_DIR}/${SAMPLE1}_sorted.bam \
  -SORT_ORDER queryname \
  -VALIDATION_STRINGENCY LENIENT \
  -MAX_RECORDS_IN_RAM 5000000 \
  -ID 1 \
  -LB CTTGTA \
  -PL ILLUMINA \
  -PU 1234 \
  -SM jurkat

rm ${OUTPUT_DIR}/${SAMPLE1}_aln.sam


### Sequence-marking duplicates
samtools fixmate -m -@ $THREADS \
  ${OUTPUT_DIR}/${SAMPLE1}_sorted.bam \
  ${OUTPUT_DIR}/${SAMPLE1}_fixmate.bam
samtools sort -@ $THREADS -m 6G \
  -o ${OUTPUT_DIR}/${SAMPLE1}_sorted.bam \
  ${OUTPUT_DIR}/${SAMPLE1}_fixmate.bam
samtools markdup -s -@ $THREADS \
  ${OUTPUT_DIR}/${SAMPLE1}_sorted.bam \
  ${OUTPUT_DIR}/${SAMPLE1}_sorted_dedup.bam
samtools index -@ $THREADS \
  ${OUTPUT_DIR}/${SAMPLE1}_sorted_dedup.bam


# ## picard mark duplicates (is this marking all of them with `PG:Z:MarkDuplicates` ?)
# java -jar /home/vbarbo/programs/picard/build/libs/picard.jar MarkDuplicates \
#   -I ${OUTPUT_DIR}/${SAMPLE1}_sorted.bam \
#   -O ${OUTPUT_DIR}/${SAMPLE1}_dedup.bam \
#   -M ${OUTPUT_DIR}/${SAMPLE1}_metrics.txt \
#   -ASSUME_SORTED true \
#   -VALIDATION_STRINGENCY LENIENT


# ## index the bam file
# java -jar /home/vbarbo/programs/picard/build/libs/picard.jar BuildBamIndex \
#   -I ${OUTPUT_DIR}/${SAMPLE1}_dedup.bam






##########################################
###### jurkat_wgs_pe_150bp ###############
###### align reads to the reference ######
##########################################


### inputs
READ1=/home/vbarbo/project_2021/paper_analysis/jurkat/data/dna_short_reads/jurkat_wgs_pe_150bp/jurkat_wgs_pe_150bp_1.fastq.gz
READ2=/home/vbarbo/project_2021/paper_analysis/jurkat/data/dna_short_reads/jurkat_wgs_pe_150bp/jurkat_wgs_pe_150bp_2.fastq.gz
SAMPLE2=jurkat150bp


### align reads to the reference
bwa mem -t $THREADS -Ma \
  -R @RG\\tID:${SAMPLE2}\\tSM:${SAMPLE2}\\tPL:ILM\\tLB:${SAMPLE2} \
  $REF \
  $READ1 \
  $READ2 \
  > ${OUTPUT_DIR}/${SAMPLE2}_aln.sam


### create read group
#java -XX:ParallelGCThreads=$THREADS -jar /home/vbarbo/programs/picard/build/libs/picard.jar AddOrReplaceReadGroups \
java -jar /home/vbarbo/programs/picard/build/libs/picard.jar AddOrReplaceReadGroups \
  -I ${OUTPUT_DIR}/${SAMPLE2}_aln.sam \
  -O ${OUTPUT_DIR}/${SAMPLE2}_sorted.bam \
  -SORT_ORDER queryname \
  -VALIDATION_STRINGENCY LENIENT \
  -ID 1 \
  -LB CTTGTA \
  -PL ILLUMINA \
  -PU 1234 \
  -SM jurkat

rm ${OUTPUT_DIR}/${SAMPLE2}_aln.sam


### Sequence-marking duplicates
samtools fixmate -m -@ $THREADS \
  ${OUTPUT_DIR}/${SAMPLE2}_sorted.bam \
  ${OUTPUT_DIR}/${SAMPLE2}_fixmate.bam
samtools sort -@ $THREADS -m 6G \
  -o ${OUTPUT_DIR}/${SAMPLE2}_sorted.bam \
  ${OUTPUT_DIR}/${SAMPLE2}_fixmate.bam
samtools markdup -s -@ $THREADS \
  ${OUTPUT_DIR}/${SAMPLE2}_sorted.bam \
  ${OUTPUT_DIR}/${SAMPLE2}_sorted_dedup.bam
samtools index -@ $THREADS \
  ${OUTPUT_DIR}/${SAMPLE2}_sorted_dedup.bam






#########################################################################################
###### merge dna_100 and dna_150 bam files to call variants with high(er) coverage ######
#########################################################################################


### inputs
SAMPLE12=jurkat_merged

DBSNP=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/Homo_sapiens_assembly38.known_indels.vcf.gz
OMNI=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/1000G_omni2.5.hg38.vcf.gz
G1000=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz
HAPMAP=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/hapmap_3.3.hg38.vcf.gz
MILLS=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz




### merge bams
samtools merge \
  -@ $THREADS \
  ${OUTPUT_DIR}/${SAMPLE12}.bam \
  ${OUTPUT_DIR}/${SAMPLE1}_sorted_dedup.bam \
  ${OUTPUT_DIR}/${SAMPLE2}_sorted_dedup.bam


### remove non-primary aligments (i.e., supplementary and secundary)
samtools view -S -b -F 2308 -@ $THREADS \
  -o ${OUTPUT_DIR}/${SAMPLE12}_primary.bam \
  ${OUTPUT_DIR}/${SAMPLE12}.bam


### index bam
samtools index -@ $THREADS \
  ${OUTPUT_DIR}/${SAMPLE12}_primary.bam







###############################################
###### call variants from the merged bam ######
###############################################


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
ref_interval_dir=/home/vbarbo/project_2021/paper_analysis/reference/genome/interval_list
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" ScatterIntervalsByNs \
  -R $REF \
  -O $ref_interval_dir/ref.interval_list
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" SplitIntervals \
  -R $REF \
  -L $ref_interval_dir/ref.interval_list \
  --scatter-count $THREADS \
  -O $ref_interval_dir/ref.scattered.interval_list

# for the next script that also use the same reference genome and the same number o threads, there is no need to create these files again.
scattered_interval_list=$ref_interval_dir/ref.scattered.interval_list


loop_num=`expr $THREADS - 1`


### calculate base recalibration
mkdir ${OUTPUT_DIR}/base_recalibration_tables
for i in `seq -f '%04g' 0 $loop_num`
do
  gatk --java-options "-Xmx4G" BaseRecalibrator \
    -R $REF \
    -I ${OUTPUT_DIR}/${SAMPLE12}_primary.bam \
    -O ${OUTPUT_DIR}/base_recalibration_tables/${SAMPLE12}_recal_data_$i.table \
    -L $scattered_interval_list/$i-scattered.interval_list \
    --known-sites $DBSNP \
    --known-sites $MILLS \
    --known-sites $G1000 &
#    --maximum-cycle-value 100000 \    # first, try without `--maximum-cycle-value 100000`
done
wait


### apply base recalibration
mkdir ${OUTPUT_DIR}/base_recalibration_bams
for i in `seq -f '%04g' 0 $loop_num`
do
  gatk --java-options "-Xmx4G" ApplyBQSR \
  -R $REF \
  -I ${OUTPUT_DIR}/${SAMPLE12}_primary.bam \
  -bqsr ${OUTPUT_DIR}/base_recalibration_tables/${SAMPLE12}_recal_data_$i.table \
  -O ${OUTPUT_DIR}/base_recalibration_bams/${SAMPLE12}_recal_$i.bam \
  -L $scattered_interval_list/$i-scattered.interval_list \
  --static-quantized-quals 10 \
  --static-quantized-quals 20 \
  --static-quantized-quals 30 &
done
wait



### Variant calling using GATK HaplotypeCaller with GVCF mode
mkdir ${OUTPUT_DIR}/intermediate_gvcfs
for i in `seq -f '%04g' 0 $loop_num`
do
  gatk --java-options "-Xmx4G" HaplotypeCaller \
    -R $REF \
    -I ${OUTPUT_DIR}/base_recalibration_bams/${SAMPLE12}_recal_$i.bam \
    -O ${OUTPUT_DIR}/intermediate_gvcfs/${SAMPLE12}_recal_$i.g.vcf \
    -L $scattered_interval_list/$i-scattered.interval_list \
    --native-pair-hmm-threads 1 \
    -ERC GVCF &
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
  -V ${OUTPUT_DIR}/intermediate_gvcfs/${SAMPLE12}_recal_$i.g.vcf \
  -L $scattered_interval_list/$i-scattered.interval_list \
  -O ${OUTPUT_DIR}/vcf_parts/${SAMPLE12}_variants_$i.vcf &
done
wait


# merge scattered phenotype vcf files
vcfFilesArg=(${OUTPUT_DIR}/vcf_parts/${SAMPLE12}_variants_*.vcf)
vcfFilesArg1=${vcfFilesArg[@]/#/-I }
gatk --java-options "-Xmx4G" GatherVcfs \
  -R $REF \
  $vcfFilesArg1 \
  -O ${OUTPUT_DIR}/${SAMPLE12}.vcf


### VARIANT QUALITY SCORE RECALIBRATION - SNPs
mkdir ${OUTPUT_DIR}/variant_recalibration
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" VariantRecalibrator \
  -V ${OUTPUT_DIR}/${SAMPLE12}.vcf \
  -O ${OUTPUT_DIR}/variant_recalibration/${SAMPLE12}_recalibrate_SNP.recal \
  -mode SNP \
  --tranches-file ${OUTPUT_DIR}/variant_recalibration/${SAMPLE12}_recalibrate_SNP.tranches \
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
  -V ${OUTPUT_DIR}/${SAMPLE12}.vcf \
  -O ${OUTPUT_DIR}/variant_recalibration/${SAMPLE12}_recalibrated_snps_raw_indels.vcf \
  --recal-file ${OUTPUT_DIR}/variant_recalibration/${SAMPLE12}_recalibrate_SNP.recal \
  --tranches-file ${OUTPUT_DIR}/variant_recalibration/${SAMPLE12}_recalibrate_SNP.tranches \
  -truth-sensitivity-filter-level 99.5 \
  --create-output-variant-index true \
  -mode SNP


### Run Variant Recalibrator – Indels
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" VariantRecalibrator \
  -V ${OUTPUT_DIR}/variant_recalibration/${SAMPLE12}_recalibrated_snps_raw_indels.vcf \
  -O ${OUTPUT_DIR}/variant_recalibration/${SAMPLE12}_recalibrate_INDEL.recal \
  -mode INDEL \
  --tranches-file ${OUTPUT_DIR}/variant_recalibration/${SAMPLE12}_recalibrate_INDEL.tranches \
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
  -V ${OUTPUT_DIR}/variant_recalibration/${SAMPLE12}_recalibrated_snps_raw_indels.vcf \
  -O ${OUTPUT_DIR}/${SAMPLE12}.recal.vcf \
  --recal-file ${OUTPUT_DIR}/variant_recalibration/${SAMPLE12}_recalibrate_INDEL.recal \
  --tranches-file ${OUTPUT_DIR}/variant_recalibration/${SAMPLE12}_recalibrate_INDEL.tranches \
  -truth-sensitivity-filter-level 99.0 \
  --create-output-variant-index true \
  -mode INDEL


### filter only PASS
bcftools view -f PASS \
  ${OUTPUT_DIR}/${SAMPLE12}.recal.vcf \
  > ${OUTPUT_DIR}/${SAMPLE12}.recal_pass.vcf


### compress and index the ground-truth VCF
bgzip -c ${OUTPUT_DIR}/${SAMPLE12}.recal_pass.vcf \
  > ${OUTPUT_DIR}/${SAMPLE12}.recal_pass.vcf.gz
rm ${OUTPUT_DIR}/${SAMPLE12}.recal_pass.vcf
bcftools index ${OUTPUT_DIR}/${SAMPLE12}.recal_pass.vcf.gz
tabix -p vcf ${OUTPUT_DIR}/${SAMPLE12}.recal_pass.vcf.gz


### delete temporary files
rm -r \
  ${OUTPUT_DIR}/base_recalibration_bams/ \
  ${OUTPUT_DIR}/base_recalibration_tables/ \
  ${OUTPUT_DIR}/intermediate_gvcfs/ \
  ${OUTPUT_DIR}/variant_recalibration/ \
  ${OUTPUT_DIR}/vcf_parts/

