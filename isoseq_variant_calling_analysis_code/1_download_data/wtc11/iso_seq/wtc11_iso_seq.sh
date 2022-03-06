# Iso-Seq reads (PacBio BAM) is downloaded from https://www.ncbi.nlm.nih.gov/sra/SRX14278671[accn]

cd /home/vbarbo/project_2021/paper_analysis/wtc11/data/iso_seq

wget https://sra-download.ncbi.nlm.nih.gov/traces/sra33/SRR/017705/SRR18130587

sam-dump --unaligned SRR18130587 > XGSUV_20200804_S64049_PL100158447-1_A01.ccs.sam

samtools view -S -b XGSUV_20200804_S64049_PL100158447-1_A01.ccs.sam > XGSUV_20200804_S64049_PL100158447-1_A01.ccs.bam

rm SRR18130587 \
  XGSUV_20200804_S64049_PL100158447-1_A01.ccs.sam



# To remove primers with limma, download primer fasta

wget https://zenodo.org/record/6332914/files/NEB_primers.fasta
