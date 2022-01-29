REF=/home/vbarbo/project_2021/paper_analysis/reference/genome/GRCh38.p13_all_chr.fasta
INPUT_BAM_DIR=/home/vbarbo/project_2021/paper_analysis/jurkat/data_manipulation
OUTPUT_DIR=/home/vbarbo/project_2021/paper_analysis/jurkat/variant_calling_from_isoseq/nanocaller
THREADS=50




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


### delete temporary files
rm -r \
  $OUTPUT_DIR/nc/intermediate_files/ \
  $OUTPUT_DIR/nc/logs/ \
  $OUTPUT_DIR/nc_sncr/intermediate_files/ \
  $OUTPUT_DIR/nc_sncr/logs/ \
  $OUTPUT_DIR/nc_sncr_fc/intermediate_files/ \
  $OUTPUT_DIR/nc_sncr_fc/logs/

