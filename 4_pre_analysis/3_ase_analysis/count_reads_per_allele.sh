#!/bin/bash




#########################
###### align reads ######
#########################

REF_STAR_INDEX=/home/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/star_indexes/star_index_sjdbOverhang99
REF_FASTA=/home/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta
GTF=/home/vbarbo/project_2021/datasets/reference/gtf_GRCh38.p13/gencode.v38.annotation.gtf
SJDBOVERHANG=99
THREADS=60
S1R1=/home/vbarbo/project_2021/datasets/wtc11/data/rna_reads/illumina/GSE175048/SRR14637256.1_1.fastq.gz
S1R2=/home/vbarbo/project_2021/datasets/wtc11/data/rna_reads/illumina/GSE175048/SRR14637256.1_2.fastq.gz
S2R1=/home/vbarbo/project_2021/datasets/wtc11/data/rna_reads/illumina/GSE175048/SRR14637257.1_1.fastq.gz
S2R2=/home/vbarbo/project_2021/datasets/wtc11/data/rna_reads/illumina/GSE175048/SRR14637257.1_2.fastq.gz
S3R1=/home/vbarbo/project_2021/datasets/wtc11/data/rna_reads/illumina/GSE175048/SRR14637258.1_1.fastq.gz
S3R2=/home/vbarbo/project_2021/datasets/wtc11/data/rna_reads/illumina/GSE175048/SRR14637258.1_2.fastq.gz
S4R1=/home/vbarbo/project_2021/datasets/wtc11/data/rna_reads/illumina/ENCSR673UKZ/ENCLB122OCH_1.fastq.gz
S4R2=/home/vbarbo/project_2021/datasets/wtc11/data/rna_reads/illumina/ENCSR673UKZ/ENCLB122OCH_2.fastq.gz
S5R1=/home/vbarbo/project_2021/datasets/wtc11/data/rna_reads/illumina/ENCSR673UKZ/ENCLB366GPZ_1.fastq.gz
S5R2=/home/vbarbo/project_2021/datasets/wtc11/data/rna_reads/illumina/ENCSR673UKZ/ENCLB366GPZ_2.fastq.gz
S6R1=/home/vbarbo/project_2021/datasets/wtc11/data/rna_reads/illumina/ENCSR673UKZ/ENCLB979NPE_1.fastq.gz
S6R2=/home/vbarbo/project_2021/datasets/wtc11/data/rna_reads/illumina/ENCSR673UKZ/ENCLB979NPE_2.fastq.gz
OUT_PREFIX=/home/vbarbo/project_2021/datasets/wtc11/ase_analysis/GSE175048_ENCSR673UKZ/star_aln/GSE175048_ENCSR673UKZ_



### star genome indexing
STAR \
  --runMode genomeGenerate \
  --genomeDir ${REF_STAR_INDEX} \
  --genomeFastaFiles ${REF_FASTA} \
  --sjdbGTFfile ${GTF} \
  --sjdbOverhang ${SJDBOVERHANG} \
  --runThreadN ${THREADS}



### this is to fix the error: "Also check ulimit -n and increase it to allow more open files"
ulimit -n 100000

### align read
STAR \
  --runMode alignReads \
  --genomeDir ${REF_STAR_INDEX} \
  --readFilesIn ${S1R1},${S2R1},${S3R1},${S4R1},${S5R1},${S6R1} \
                ${S1R2},${S2R2},${S3R2},${S4R2},${S5R2},${S6R2} \
  --readFilesCommand zcat \
  --outFileNamePrefix ${OUT_PREFIX} \
  --outSAMmapqUnique 60 \
  --runThreadN ${THREADS} \
  --sjdbOverhang ${SJDBOVERHANG} \
  --bamRemoveDuplicatesType - \
  --twopassMode Basic \
  --outSAMtype BAM SortedByCoordinate \
  --sjdbGTFfile ${GTF}


### remove non-primoary alignments: unmapped (4), secondary (256), supplementary (2048)
### don't remove duplicate (1024)
samtools view \
  -S -b \
  -F 2308 \
  -@ $THREADS \
  -o ${OUT_PREFIX}Aligned.sortedByCoord.out_primary.bam \
  ${OUT_PREFIX}Aligned.sortedByCoord.out.bam




######################
### add read group ###
######################

INPUT_BAM=/home/vbarbo/project_2021/datasets/wtc11/ase_analysis/GSE175048_ENCSR673UKZ/star_aln/GSE175048_ENCSR673UKZ_Aligned.sortedByCoord.out_primary.bam
output_bam=${INPUT_BAM%.bam}_readGroupAddedSameOrder.bam

#java -jar /home/vbarbo/programs/picard/build/libs/picard.jar AddOrReplaceReadGroups \
java -XX:ParallelGCThreads=$THREADS -jar /home/vbarbo/programs/picard/build/libs/picard.jar AddOrReplaceReadGroups \
  -I ${INPUT_BAM} \
  -O ${output_bam} \
  -VALIDATION_STRINGENCY LENIENT \
  -MAX_RECORDS_IN_RAM 5000000 \
  -ID 1 \
  -LB lib1 \
  -PL ILLUMINA \
  -PU unit1 \
  -SM wtc11
#  -SORT_ORDER coordinate
#  -SORT_ORDER queryname


### index bam
samtools index \
  -@ ${THREADS} \
  ${output_bam}



########################################################################################################
###### Calculate read counts per allele for allele-specific expression analysis of RNAseq data #########
########################################################################################################
# https://gatk.broadinstitute.org/hc/en-us/articles/360037054312-ASEReadCounter
# compatible with Mamba, a downstream tool developed for allele-specific expression analysis
# This tool will only process biallelic het SNP sites


REF_FASTA=/home/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta
SITES_VCF=/home/vbarbo/project_2021/datasets/wtc11/truth/allen_institute_wtc11_shortreads_wgs/3546dc62_AH77TTBBXX_DS-229105_GCCAAT_recalibrated_subsetChromosomes_pass.vcf.gz
INPUT_BAM=/home/vbarbo/project_2021/datasets/wtc11/ase_analysis/GSE175048_ENCSR673UKZ/star_aln/GSE175048_ENCSR673UKZ_Aligned.sortedByCoord.out_primary_readGroupAddedSameOrder.bam
OUTPUT_TABLE_DIR=/home/vbarbo/project_2021/datasets/wtc11/ase_analysis/GSE175048_ENCSR673UKZ/ase_read_count_tables



### keep only snps
vcftools --gzvcf ${SITES_VCF} \
  --out ${SITES_VCF%.vcf.gz}_snps \
  --remove-indels \
  --recode \
  --recode-INFO-all


### compressing and indexing
sites_vc_snp=${SITES_VCF%.vcf.gz}_snps.recode.vcf.gz
bgzip -c \
  ${SITES_VCF%.vcf.gz}_snps.recode.vcf \
  ${sites_vc_snp}
bcftools index ${sites_vc_snp}
tabix -p vcf ${sites_vc_snp}
rm ${SITES_VCF%.vcf.gz}_snps.recode.vcf





SCATTERED_INTERVAL_LIST=/home/vbarbo/project_2021/datasets/gloria_data/analysis/my_ground_truth_jurkat_wgs_pe_100bp/jurkat100bp_scattered.interval_list
THREADS=30
loop_num=`expr $THREADS - 1`

mkdir ${OUTPUT_TABLE_DIR}

for i in `seq -f '%04g' 0 $loop_num`
do
  gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" ASEReadCounter \
    -R ${REF_FASTA} \
    -I ${INPUT_BAM} \
    -V ${sites_vc_snp} \
    -O ${OUTPUT_TABLE_DIR}/ase_read_count_$i.table \
    -L $SCATTERED_INTERVAL_LIST/$i-scattered.interval_list &
done
wait



### concatenate the table files to a single table file
table_files=(`ls ${OUTPUT_TABLE_DIR}/ase_read_count_0*.table`)

cat ${table_files[0]} > ${OUTPUT_TABLE_DIR}/ase_read_count_all.table

for x in `echo ${table_files[@]:1}`
do
  sed '1 d' `echo ${x}` >> ${OUTPUT_TABLE_DIR}/ase_read_count_all.table
done
wait


