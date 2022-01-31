# data downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175048


cd /home/vbarbo/project_2021/paper_analysis/wtc11/data/rna_short_reads/GSE175048

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-25/SRR14637256/SRR14637256.1
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-25/SRR14637257/SRR14637257.1
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-25/SRR14637258/SRR14637258.1

fastq-dump --split-files SRR14637256.1
fastq-dump --split-files SRR14637257.1
fastq-dump --split-files SRR14637258.1

bgzip -c SRR14637256.1_1.fastq > SRR14637256.1_1.fastq.gz
bgzip -c SRR14637256.1_2.fastq > SRR14637256.1_2.fastq.gz

bgzip -c SRR14637257.1_1.fastq > SRR14637257.1_1.fastq.gz
bgzip -c SRR14637257.1_2.fastq > SRR14637257.1_2.fastq.gz

bgzip -c SRR14637258.1_1.fastq > SRR14637258.1_1.fastq.gz
bgzip -c SRR14637258.1_2.fastq > SRR14637258.1_2.fastq.gz

rm SRR14637256.1 SRR14637257.1 SRR14637258.1 \
  SRR14637256.1_1.fastq SRR14637256.1_2.fastq \
  SRR14637257.1_1.fastq SRR14637257.1_2.fastq \
  SRR14637258.1_1.fastq SRR14637258.1_2.fastq

