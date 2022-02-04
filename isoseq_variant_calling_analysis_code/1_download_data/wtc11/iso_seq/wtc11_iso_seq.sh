# Download the Iso-Seq CCS reads for WTC-11
wget https://sheynkman-lab-lifebit.s3.eu-west-1.amazonaws.com/Sample-Data/WTC11/XGSUV_20200804_S64049_PL100158447-1_A01.ccs.bam \
  -P /home/vbarbo/project_2021/paper_analysis/wtc11/data/iso_seq/

# To remove primers, we used the following FASTA files of primer sequences.
wget https://sheynkman-lab-lifebit.s3.eu-west-1.amazonaws.com/Reference/NEB_primers.fasta \
  -P /home/vbarbo/project_2021/paper_analysis/wtc11/data/iso_seq/

