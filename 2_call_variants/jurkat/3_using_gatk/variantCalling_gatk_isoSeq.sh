#!/bin/bash


# wtc-11
# iso-seq data
# read group added
# duplicates not marked/removed
# only primary alignments kept
# using sncr
# without flagcorrection



### inputs
REF=/home/vbarbo/project_2021/paper_analysis/reference/genome/GRCh38.p13_all_chr.fasta
INPUT_BAM=/home/vbarbo/project_2021/paper_analysis/jurkat/data_manipulation/aln_sncr.bam
OUTPUT_DIR=/home/vbarbo/project_2021/paper_analysis/jurkat/variant_calling_from_isoseq/gatk
SAMPLE=isoSeq_jurkat
THREADS=30

DBSNP=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/Homo_sapiens_assembly38.known_indels.vcf.gz
OMNI=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/1000G_omni2.5.hg38.vcf.gz
G1000=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz
HAPMAP=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/hapmap_3.3.hg38.vcf.gz
MILLS=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz



### to use 30 cores, intervals for the reference genome were already created in 
### /home/vbarbo/project_2021/projects/lrRNA-seq_variant_calling/1_generate_ground_truth/jurkat/generateGroundTruth_shortReads_gatk_jurkat.sh
SCATTERED_INTERVAL_LIST=/home/vbarbo/project_2021/datasets/gloria_data/analysis/my_ground_truth_jurkat_wgs_pe_100bp/jurkat100bp_scattered.interval_list
THREADS=30
loop_num=`expr $THREADS - 1`




### Add or replace read groups
# to read about read group: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4243306/
# topic: Preparing the appropriate read group information
java -XX:ParallelGCThreads=$THREADS -jar /home/vbarbo/programs/picard/build/libs/picard.jar AddOrReplaceReadGroups \
  -I $INPUT_BAM \
  -O $OUTPUT_DIR/aln_sncr_rg.bam \
  -SO coordinate \
  -RGID group1 \
  -RGSM sample1 \
  -RGPL pacbio \
  -RGLB lib1 \
  -RGPU unit1

samtools index -@ $THREADS $OUTPUT_DIR/aln_sncr_rg.bam



### calculate base recalibration
mkdir ${OUTPUT_DIR}/base_recalibration_tables
for i in `seq -f '%04g' 0 $loop_num`
do
  gatk --java-options "-Xmx4G" BaseRecalibrator \
    -R $REF \
    -I $OUTPUT_DIR/aln_sncr_rg.bam \
    -O ${OUTPUT_DIR}/base_recalibration_tables/${SAMPLE}_recal_data_$i.table \
    -L $SCATTERED_INTERVAL_LIST/$i-scattered.interval_list \
    --known-sites $DBSNP \
    --known-sites $MILLS \
    --known-sites $G1000 \
    --maximum-cycle-value 100000 &
done
wait



### apply base recalibration
mkdir ${OUTPUT_DIR}/base_recalibration_bams
for i in `seq -f '%04g' 0 $loop_num`
do
  gatk --java-options "-Xmx4G" ApplyBQSR \
  -R $REF \
  -I $OUTPUT_DIR/aln_sncr_rg.bam \
  -bqsr ${OUTPUT_DIR}/base_recalibration_tables/${SAMPLE}_recal_data_$i.table \
  -O ${OUTPUT_DIR}/base_recalibration_bams/${SAMPLE}_recal_$i.bam \
  -L $SCATTERED_INTERVAL_LIST/$i-scattered.interval_list \
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
    -I ${OUTPUT_DIR}/base_recalibration_bams/${SAMPLE}_recal_$i.bam \
    -O ${OUTPUT_DIR}/intermediate_gvcfs/${SAMPLE}_recal_$i.g.vcf \
    -L $SCATTERED_INTERVAL_LIST/$i-scattered.interval_list \
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
  -V ${OUTPUT_DIR}/intermediate_gvcfs/${SAMPLE}_recal_$i.g.vcf \
  -L $SCATTERED_INTERVAL_LIST/$i-scattered.interval_list \
  -O ${OUTPUT_DIR}/vcf_parts/${SAMPLE}_variants_$i.vcf &
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



### filter only PASS
bcftools view -f PASS \
  ${OUTPUT_DIR}/${SAMPLE}.recal.vcf \
  > ${OUTPUT_DIR}/${SAMPLE}.recal_pass.vcf



### sort, compress and index
bcftools sort \
  -O z \
  -o ${OUTPUT_DIR}/${SAMPLE}.recal_pass.vcf.gz \
  ${OUTPUT_DIR}/${SAMPLE}.recal_pass.vcf

bcftools index ${OUTPUT_DIR}/${SAMPLE}.recal_pass.vcf.gz
tabix -p vcf ${OUTPUT_DIR}/${SAMPLE}.recal_pass.vcf.gz


