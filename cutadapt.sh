#!/bin/bash
#SBATCH --job-name=cutadapt
#SBATCH --mail-user=surendra.neupane@moffitt.org
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --time=24:00:00
date;hostname;pwd

# Load the cutadapt module
module load cutadapt/2.10-GCCcore-9.3.0-Python-3.8.2

# Define the samples array with the updated sample names
samples=(
    "S94_CP3_D14_SFEM__1312024_S94_ME_L001"
    "S93_CP3_D7_SFEM__1312024_S93_ME_L001"
    "S92_CD3_D12_SFEM__1312024_S92_ME_L001"
    "S91_CD3_D7_SFEM__1312024_S91_ME_L001"
    "S90_CD2_D7_SFEM__1312024_S4S90_ME_L001"
    "S89_CP2_D7_SFEM__1312024_S89_ME_L001"
    "S88_CP1_D14_SFEM__1312024_S88_ME_L001"
    "S87_CP1_D7_SFEM__1312024_S87S3_ME_L001"
)

# Directory paths (adjust these paths as necessary)
input_dir="/share/lab_padron/CICPT_4628_AgilentPanel-February-2024/merged"
output_dir="/share/lab_padron/CICPT_4628_AgilentPanel-February-2024/merged"

# Loop through each sample and process it with cutadapt
for sample in "${samples[@]}"; do
    INPUT_R1="${input_dir}/${sample}_R1_001.fastq.gz"
    INPUT_R2="${input_dir}/${sample}_R2_001.fastq.gz"
    OUTPUT_R1="${output_dir}/${sample}_R1_trimmed.fastq.gz"
    OUTPUT_R2="${output_dir}/${sample}_R2_trimmed.fastq.gz"
    
    # Check if the input files exist, then run cutadapt
    if [[ -e $INPUT_R1 && -e $INPUT_R2 ]]; then
        cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o $OUTPUT_R1 -p $OUTPUT_R2 -m 20 $INPUT_R1 $INPUT_R2
    else
        echo "Files for $sample are missing!"
    fi
done
