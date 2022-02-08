## Download FLNC Iso-Seq sequences, PacBio BAM file, from Jurkat cells.
## The file is public available in https://doi.org/10.5281/zenodo.5920920.
## In the documentation, they state that, since the file was too big to be uploaded as a single file,
##   it was split between smaller files. Hence, to donwload the sequences, we need to download the parts 
##   of the BAM, along with the header, and finaly, concatenate all parts.



OUTPUT_DIR=/home/vbarbo/project_2021/paper_analysis/jurkat/data/iso_seq


cd ${OUTPUT_DIR}

### download the BAM parts
wget https://zenodo.org/record/5920920/files/LRPG-Manuscript-Results-results-results-jurkat-jurkat.flnc.chunk.aa.bam
wget https://zenodo.org/record/5920920/files/LRPG-Manuscript-Results-results-results-jurkat-jurkat.flnc.chunk.ab.bam
wget https://zenodo.org/record/5920920/files/LRPG-Manuscript-Results-results-results-jurkat-jurkat.flnc.chunk.ac.bam
wget https://zenodo.org/record/5920920/files/LRPG-Manuscript-Results-results-results-jurkat-jurkat.flnc.chunk.ad.bam
wget https://zenodo.org/record/5920920/files/LRPG-Manuscript-Results-results-results-jurkat-jurkat.flnc.chunk.ae.bam
wget https://zenodo.org/record/5920920/files/LRPG-Manuscript-Results-results-results-jurkat-jurkat.flnc.chunk.af.bam
wget https://zenodo.org/record/5920920/files/LRPG-Manuscript-Results-results-results-jurkat-jurkat.flnc.chunk.ag.bam
wget https://zenodo.org/record/5920920/files/LRPG-Manuscript-Results-results-results-jurkat-jurkat.flnc.chunk.ah.bam
wget https://zenodo.org/record/5920920/files/LRPG-Manuscript-Results-results-results-jurkat-jurkat.flnc.chunk.ai.bam
wget https://zenodo.org/record/5920920/files/LRPG-Manuscript-Results-results-results-jurkat-jurkat.flnc.chunk.aj.bam
wget https://zenodo.org/record/5920920/files/LRPG-Manuscript-Results-results-results-jurkat-jurkat.flnc.chunk.ak.bam
wget https://zenodo.org/record/5920920/files/LRPG-Manuscript-Results-results-results-jurkat-jurkat.flnc.chunk.al.bam
wget https://zenodo.org/record/5920920/files/LRPG-Manuscript-Results-results-results-jurkat-jurkat.flnc.chunk.am.bam

### download header
wget https://zenodo.org/record/5987905/files/LRPG-Manuscript-Results-results-results-jurkat-isoseq3-companion-files.tar.gz

tar -zxf \
  LRPG-Manuscript-Results-results-results-jurkat-isoseq3-companion-files.tar.gz \
  LRPG-Manuscript-Results/results/results/jurkat/isoseq3/jurkat.flnc.header.sam
header_sam=LRPG-Manuscript-Results/results/results/jurkat/isoseq3/jurkat.flnc.header.sam
mv ${header_sam} \
  jurkat.flnc.sam
rmdir -p $(dirname ${header_sam})

### convert the BAM parts to SAMs and join them in a single SAM
for x in $OUTPUT_DIR/LRPG-Manuscript-Results-results-results-jurkat-jurkat.flnc.chunk.a*.bam
do
  samtools view ${x} >> jurkat.flnc.sam
done

### convert SAM back to BAM
samtools view -bS \
  jurkat.flnc.sam \
  > jurkat.flnc.bam

rm jurkat.flnc.sam

