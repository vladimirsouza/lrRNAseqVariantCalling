ORI_BAM=~/project_2021/paper_analysis/wtc11/data_manipulation/aln_s_small.bam
SNCR_BAM=~/project_2021/paper_analysis/wtc11/data_manipulation/aln_sncr_small.bam
OUTPUT_DIR=~/project_2021/paper_analysis/wtc11/data_manipulation
THREADS=30



##################
ORI_BAM=~/project_2021/paper_analysis/wtc11/data_manipulation/aln_s.bam
SNCR_BAM=~/project_2021/paper_analysis/wtc11/data_manipulation/aln_sncr.bam
OUTPUT_DIR=~/project_2021/paper_analysis/wtc11/data_manipulation
THREADS=30
##################


### a temp directory to save files and delete at the end
temp_dir=$(mktemp -d ${OUTPUT_DIR}/flagCorrection_temp_dir_XXXXXXXXXXXXXXXX)


### work with pieces of the SNCR BAM to save memory
### spit BAM by chromosome
bamtools split \
  -in $SNCR_BAM \
  -stub ${temp_dir}/split_bam \
  -reference &
# real    114m18.043s
# user    111m14.526s
# sys     2m43.724s


### create an associative array with FLAGs (values) and QNAMEs (key)
# create file with QNAMEs and FLAGs of the original BAM
samtools view ${ORI_BAM} | cut --output-delimiter="=" -f 1,2 > ${temp_dir}/qname_flag.txt
# associative array from the file
declare -A -g ori_flags
readarray -t lines < ${temp_dir}/qname_flag.txt
for line in "${lines[@]}"; do
  ori_flags[${line%%=*}]=${line#*=}
done
# real    91m31.269s
# user    91m14.280s
# sys     0m15.651s


wait


### replace FLAGs of the SNCR BAM with FLAGs from the original BAM
split_sncr_bams=( ${temp_dir}/split_bam*.bam )
for i in "${!split_sncr_bams[@]}"; do
  split_sncr_sams[$i]=${split_sncr_bams[$i]%bam}sam
  bam_header[$i]=${split_sncr_bams[$i]%.bam}_header.txt
  split_sncr_qnames[$i]=${split_sncr_bams[$i]%.bam}_qnames.txt
  corrected_flags_file[$i]=${split_sncr_bams[$i]%.bam}_corrected_flags.txt
  split_sncr_sams_correct[$i]=${split_sncr_bams[$i]%.bam}_correct.sam
  split_sncr_sams_correct_header[$i]=${split_sncr_bams[$i]%bam}_correct_header.sam
  split_sncr_sams_correct_header_bam[$i]=${split_sncr_bams[$i]%bam}_correct_header.bam
done

for i in "${!split_sncr_bams[@]}"; do
  ### bam to sam
  samtools view -o ${split_sncr_sams[$i]} ${split_sncr_bams[$i]}
  samtools view -H ${split_sncr_bams[$i]} > ${bam_header[$i]}
  rm ${split_sncr_bams[$i]}
  
  ### get array with corrected flags for sncr sam
  cut -f 1 ${split_sncr_sams[$i]} > ${split_sncr_qnames[$i]}
  # how could i use multiple cores here?
  # each loop needs to create arrays with different names to not subwrite each other
  # or could i avoid arrays?
  # could i use a multidementional array: columns refer to sam files, rows refer to flags
  readarray -t qnames < ${split_sncr_qnames[$i]}
  for j in "${!qnames[@]}"; do
    t=${qnames[$j]}
    corrected_flags[$j]=${ori_flags[$t]}
  done
  ### save array to a file
  printf "%s\n" "${corrected_flags[@]}" > ${corrected_flags_file[$i]}
  
  ### replace flag column of the sam file
  awk 'FNR==NR{a[NR]=$1;next}{$2=a[FNR]}1' OFS='\t' \
    ${corrected_flags_file[$i]} \
    ${split_sncr_sams[$i]} \
    > ${split_sncr_sams_correct[$i]}
  rm ${split_sncr_sams[$i]}
  
  ### add headers to the sam
  cat ${bam_header[$i]} | cat - ${split_sncr_sams_correct[$i]} > ${split_sncr_sams_correct_header[$i]}
  rm ${split_sncr_sams_correct[$i]}
  ### convert sam to bam, and delete sam
  samtools view -S -b ${split_sncr_sams_correct_header[$i]} > ${split_sncr_sams_correct_header_bam[$i]}
  rm ${split_sncr_sams_correct_header[$i]}
  rm ${bam_header[$i]} ${split_sncr_qnames[$i]} ${corrected_flags_file[$i]}
done





### merge the corrected split bams
split_sncr_bams=( ${temp_dir}/split_bam*.bam )
#########samtools merge merged.bam *.bam   # to do




