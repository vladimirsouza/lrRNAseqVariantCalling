
### NanoCaller + SNCR + flagCorrection
REF=/home/Shared_disaster_recovery/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta
BAM=/home/Shared_disaster_recovery/vbarbo/project_2021/datasets/wtc11/manipulate_data/aln_sncr_fc.bam
OUTPUT_DIR=/home/vbarbo/project_2021/variant_calling_isoseq/wtc11/nanocaller/yes_sncr_yes_fc
THREADS=20

source activate NanoCaller

python /home/vbarbo/Programs/nanocaller_conda/NanoCaller/scripts/NanoCaller_WGS.py \
  -bam ${BAM} \
  -ref ${REF} \
  -o $OUTPUT_DIR \
  -cpu $THREADS \
  -p ccs





### NanoCaller (alone)
REF=/home/Shared_disaster_recovery/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta
BAM=/home/Shared_disaster_recovery/vbarbo/project_2021/datasets/wtc11/manipulate_data/aln_s.bam
OUTPUT_DIR=/home/vbarbo/project_2021/variant_calling_isoseq/wtc11/nanocaller/no_sncr
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
BAM=/home/Shared_disaster_recovery/vbarbo/project_2021/datasets/wtc11/manipulate_data/aln_sncr.bam
OUTPUT_DIR=/home/vbarbo/project_2021/variant_calling_isoseq/wtc11/nanocaller/yes_sncr_no_fc
THREADS=10

source activate NanoCaller

python /home/vbarbo/Programs/nanocaller_conda/NanoCaller/scripts/NanoCaller_WGS.py \
  -bam ${BAM} \
  -ref ${REF} \
  -o $OUTPUT_DIR \
  -cpu $THREADS \
  -p ccs

