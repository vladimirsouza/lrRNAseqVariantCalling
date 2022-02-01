### get the short-read fatsq files
### there are two datasets: jurkat_wgs_pe_100bp and jurkat_wgs_pe_150bp

### We download these datasets because we want to generate a high-quality diploid VCF file to use as the ground truth for Jurkat cells.



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

