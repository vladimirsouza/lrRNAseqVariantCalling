REF=/home/vbarbo/project_2021/paper_analysis/reference/genome/GRCh38.p13_all_chr.fasta
INPUT_BAM_DIR=/home/vbarbo/project_2021/paper_analysis/wtc11/data_manipulation
OUTPUT_DIR=/home/vbarbo/project_2021/paper_analysis/wtc11/variant_calling_from_isoseq/nanocaller
THREADS=30



# ### NanoCaller (alone)
# singularity exec -e --pwd /app \
#   /home/vbarbo/programs/NanoCaller/nanocaller-1.0.1.simg \
#   python NanoCaller.py \
#   -bam $INPUT_BAM_DIR/aln.bam \
#   -p ccs \
#   -o $OUTPUT_DIR/nc \
#   -ref $REF \
#   -cpu $THREADS


conda deactivate
conda activate NanoCaller


### NanoCaller (alone)
python /home/vbarbo/programs/NanoCaller/scripts/NanoCaller_WGS.py \
  -bam $INPUT_BAM_DIR/aln_s.bam \
  -ref $REF \
  -o $OUTPUT_DIR/nc \
  -cpu $THREADS \
  -p ccs

### NanoCaller + SNCR
python /home/vbarbo/programs/NanoCaller/scripts/NanoCaller_WGS.py \
  -bam $INPUT_BAM_DIR/aln_sncr.bam \
  -ref $REF \
  -o $OUTPUT_DIR/nc_sncr \
  -cpu $THREADS \
  -p ccs

### NanoCaller + SNCR + flagCorrection
python /home/vbarbo/programs/NanoCaller/scripts/NanoCaller_WGS.py \
  -bam $INPUT_BAM_DIR/aln_sncr_fc.bam \
  -ref $REF \
  -o $OUTPUT_DIR/nc_sncr_fc \
  -cpu $THREADS \
  -p ccs


