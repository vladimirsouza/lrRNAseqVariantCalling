#!/bin/bash



############################
###### install clair3 ######
############################


conda deactivate

# create and activate an environment named clair3
conda create -n clair3 python=3.6.10 -y
source activate clair3

# install pypy and packages in the environemnt
conda install -c conda-forge pypy3.6 -y
pypy3 -m ensurepip
pypy3 -m pip install mpmath==1.2.1

# install python packages in environment
pip3 install tensorflow==2.2.0
pip3 install tensorflow-addons==0.11.2 tables==3.6.1
conda install -c anaconda pigz==2.4 -y
conda install -c conda-forge parallel=20191122 zstd=1.4.4 -y
conda install -c conda-forge -c bioconda samtools=1.10 -y
conda install -c conda-forge -c bioconda whatshap=1.0 -y

# clone Clair3
git clone https://github.com/HKU-BAL/Clair3.git
cd /home/vbarbo/programs/Clair3

# download pre-trained models
mkdir models
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz 
tar -zxvf clair3_models.tar.gz -C ./models



###################################
###### run clair3 on iso-seq ######
###################################


conda deactivate
conda activate clair3


REF=/home/vbarbo/project_2021/paper_analysis/reference/genome/GRCh38.p13_all_chr.fasta
INPUT_BAM_DIR=/home/vbarbo/project_2021/paper_analysis/wtc11/data_manipulation
OUTPUT_DIR=/home/vbarbo/project_2021/paper_analysis/wtc11/variant_calling_from_isoseq/clair3
THREADS=20
PATH_TO_REPO=/home/vbarbo/project_2021/projects/lrRNAseqVariantCalling




### Clair3 (alone)
# call variants
/home/vbarbo/programs/Clair3/run_clair3.sh \
  --bam_fn=$INPUT_BAM_DIR/aln_s.bam \
  --ref_fn=$REF \
  --threads=$THREADS \
  --platform="hifi" \
  --model_path="/home/vbarbo/programs/Clair3/models/hifi" \
  --output=$OUTPUT_DIR/c3
# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/c3/pileup.vcf.gz \
  > $OUTPUT_DIR/c3/pileup_pass.vcf
# compress and index
bgzip -c $OUTPUT_DIR/c3/pileup_pass.vcf \
  > $OUTPUT_DIR/c3/pileup_pass.vcf.gz
bcftools index $OUTPUT_DIR/c3/pileup_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/c3/pileup_pass.vcf.gz
rm $OUTPUT_DIR/c3/pileup_pass.vcf



### SNCR + Clair3
/home/vbarbo/programs/Clair3/run_clair3.sh \
  --bam_fn=$INPUT_BAM_DIR/aln_sncr.bam \
  --ref_fn=$REF \
  --threads=$THREADS \
  --platform="hifi" \
  --model_path="/home/vbarbo/programs/Clair3/models/hifi" \
  --output=$OUTPUT_DIR/c3_sncr
# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/c3_sncr/pileup.vcf.gz \
  > $OUTPUT_DIR/c3_sncr/pileup_pass.vcf
# compress and index
bgzip -c $OUTPUT_DIR/c3_sncr/pileup_pass.vcf \
  > $OUTPUT_DIR/c3_sncr/pileup_pass.vcf.gz
bcftools index $OUTPUT_DIR/c3_sncr/pileup_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/c3_sncr/pileup_pass.vcf.gz
rm $OUTPUT_DIR/c3_sncr/pileup_pass.vcf



### SNCR + FC + Clair3
/home/vbarbo/programs/Clair3/run_clair3.sh \
  --bam_fn=$INPUT_BAM_DIR/aln_sncr_fc.bam \
  --ref_fn=$REF \
  --threads=$THREADS \
  --platform="hifi" \
  --model_path="/home/vbarbo/programs/Clair3/models/hifi" \
  --output=$OUTPUT_DIR/c3_sncr_fc
# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/c3_sncr_fc/pileup.vcf.gz \
  > $OUTPUT_DIR/c3_sncr_fc/pileup_pass.vcf
# compress and index
bgzip -c $OUTPUT_DIR/c3_sncr_fc/pileup_pass.vcf \
  > $OUTPUT_DIR/c3_sncr_fc/pileup_pass.vcf.gz
bcftools index $OUTPUT_DIR/c3_sncr_fc/pileup_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/c3_sncr_fc/pileup_pass.vcf.gz
rm $OUTPUT_DIR/c3_sncr_fc/pileup_pass.vcf



### Clair3-mix:
### SNPs from Clair3 (alone), indels from SNCR + FC+ Clair3

# snps
vcftools --gzvcf $OUTPUT_DIR/c3/pileup_pass.vcf.gz \
  --out $OUTPUT_DIR/mix/pileup_pass_snp \
  --remove-indels --recode --recode-INFO-all
bgzip $OUTPUT_DIR/mix/pileup_pass_snp.recode.vcf
tabix -p vcf $OUTPUT_DIR/mix/pileup_pass_snp.recode.vcf.gz

# indels
vcftools --gzvcf $OUTPUT_DIR/c3_sncr_fc/pileup_pass.vcf.gz \
  --out $OUTPUT_DIR/mix/pileup_pass_indel \
  --keep-only-indels --recode --recode-INFO-all
bgzip $OUTPUT_DIR/mix/pileup_pass_indel.recode.vcf
tabix -p vcf $OUTPUT_DIR/mix/pileup_pass_indel.recode.vcf.gz

# concatenate vcf files
bcftools concat \
  $OUTPUT_DIR/mix/pileup_pass_snp.recode.vcf.gz \
  $OUTPUT_DIR/mix/pileup_pass_indel.recode.vcf.gz \
  -o $OUTPUT_DIR/mix/pileup_pass_mix.recode.vcf.gz \
  -O z -D -a

# in case the concatenate vcf ends up with two different variants at a same site,
# keeping the one that shows the highest QUAL value.
Rscript $PATH_TO_REPO/removeRepeatedLowerQualSites.r \
  $OUTPUT_DIR/mix/pileup_pass_mix.recode.vcf.gz \
  $OUTPUT_DIR/mix/pileup_pass_mix_norep.recode.vcf.gz



### delete temporary files
rm -r \
  $OUTPUT_DIR/c3/log/ \
  $OUTPUT_DIR/c3/tmp/ \
  $OUTPUT_DIR/c3_sncr/log/ \
  $OUTPUT_DIR/c3_sncr/tmp/ \
  $OUTPUT_DIR/c3_sncr_fc/log/ \
  $OUTPUT_DIR/c3_sncr_fc/tmp/

