
### NanoCaller + SNCR + flagCorrection
REF=/home/Shared_disaster_recovery/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta
BAM=/home/Shared_disaster_recovery/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/noMarkDuplicate/aln_split_flagCorrection.bam
OUTPUT_DIR=/home/vbarbo/project_2021/variant_calling_isoseq/jurkat/nanocaller/yes_sncr_yes_fc
THREADS=10

source activate NanoCaller

python /home/vbarbo/Programs/nanocaller_conda/NanoCaller/scripts/NanoCaller_WGS.py \
  -bam ${BAM} \
  -ref ${REF} \
  -o $OUTPUT_DIR \
  -cpu $THREADS \
  -p ccs





### NanoCaller (alone)
REF=/home/Shared_disaster_recovery/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta
BAM=/home/Shared_disaster_recovery/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/aln.bam
OUTPUT_DIR=/home/vbarbo/project_2021/variant_calling_isoseq/jurkat/nanocaller/no_sncr
THREADS=10

source activate NanoCaller

python /home/vbarbo/Programs/nanocaller_conda/NanoCaller/scripts/NanoCaller_WGS.py \
  -bam ${BAM} \
  -ref ${REF} \
  -o $OUTPUT_DIR \
  -cpu $THREADS \
  -p ccs





### NanoCaller + SNCR
REF=/home/Shared_disaster_recovery/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta
BAM=/home/Shared_disaster_recovery/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/noMarkDuplicate/aln_split.bam
OUTPUT_DIR=/home/vbarbo/project_2021/variant_calling_isoseq/jurkat/nanocaller/yes_sncr_no_fc
THREADS=10

source activate NanoCaller

python /home/vbarbo/Programs/nanocaller_conda/NanoCaller/scripts/NanoCaller_WGS.py \
  -bam ${BAM} \
  -ref ${REF} \
  -o $OUTPUT_DIR \
  -cpu $THREADS \
  -p ccs

