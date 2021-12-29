#!/bin/bash



####### read groups added (needed for gatk)
####### no duplicates removal


### inputs
REF=/home/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta
OUTPUT_DIR=/home/vbarbo/project_2021/datasets/gloria_data/analysis/gatk_calls_isoSeq/noMarkDuplicate
SAMPLE=isoSeq_jurkat
THREADS=30

DBSNP=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/Homo_sapiens_assembly38.known_indels.vcf.gz
OMNI=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/1000G_omni2.5.hg38.vcf.gz
G1000=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz
HAPMAP=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/hapmap_3.3.hg38.vcf.gz
MILLS=/home/vbarbo/project_2021/datasets/reference/vcf_human_ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

INPUT_BAM=/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/aln.bam


# use the same intervals already created for this reference and 30 cores
# this folder was generates by script /home/vbarbo/project_2021/scripts/repository/jurkat/generate_ground_truth/generateGroundTruth_shortReads_gatk_jurkat.sh
# look for the gatk functions ScatterIntervalsByNs and SplitIntervals
SCATTERED_INTERVAL_LIST=/home/vbarbo/project_2021/datasets/gloria_data/analysis/my_ground_truth_jurkat_wgs_pe_100bp/jurkat100bp_scattered.interval_list
THREADS=30
loop_num=`expr $THREADS - 1`




### Add or replace read groups
java -XX:ParallelGCThreads=$THREADS -jar /home/vbarbo/programs/picard/build/libs/picard.jar AddOrReplaceReadGroups \
  -I $INPUT_BAM \
  -O $OUTPUT_DIR/aln_rg.bam \
  -SO coordinate \
  -RGID 4 \
  -RGLB lib1 \
  -RGPL pacbio \
  -RGPU unit1 \
  -RGSM sample1


### Split N
### there is no need to reassign the mapping quality:
### (https://lh3.github.io/minimap2/minimap2.html):
### "Mapping quality (0-255 with 255 for missing)"
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" SplitNCigarReads \
  -R $REF \
  -I $OUTPUT_DIR/aln_rg.bam \
  -O $OUTPUT_DIR/aln_rg_split.bam

samtools index -@ $THREADS $OUTPUT_DIR/aln_rg_split.bam


### calculate base recalibration
mkdir ${OUTPUT_DIR}/base_recalibration_tables
for i in `seq -f '%04g' 0 $loop_num`
do
  gatk --java-options "-Xmx4G" BaseRecalibrator \
    --maximum-cycle-value 100000 \
    -R $REF \
    -I $OUTPUT_DIR/aln_rg_split.bam \
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
  -I $OUTPUT_DIR/aln_rg_split.bam \
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




