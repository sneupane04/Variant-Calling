#!/bin/bash
#SBATCH --job-name=align_reads
#SBATCH --mail-user=surendra.neupane@moffitt.org
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --time=24:00:00
date; hostname; pwd

# Load the necessary modules
module load BWA/0.7.17-GCC-10.2.0 
module load SAMtools/1.16.1-GCC-11.3.0

# Define the human genome reference path
HUMAN_GENOME="/share/lab_padron/CICPT-4562-Padron_Traci-QIAseq-Dec-01-2023_Data/fastq/fastqtraci/GRCh38.p14.genome.fa"

# Directory for input and output
input_dir="/share/lab_padron/CICPT_4628_AgilentPanel-February-2024/merged"
output_dir="/share/lab_padron/CICPT_4628_AgilentPanel-February-2024/merged"

# Sample names derived from provided filenames
samples=(
    "S96_Enzo-Positive-Control--20ng_Enzo-Positive-Control--20ng_1312024_S6S96_ME"
    "S87_CP1_D7_SFEM__1312024_S87S3_ME"
    "S90_CD2_D7_SFEM__1312024_S4S90_ME"
    "S88_CP1_D14_SFEM__1312024_S88_ME_L001"
    "S89_CP2_D7_SFEM__1312024_S89_ME_L001"
    "S91_CD3_D7_SFEM__1312024_S91_ME_L001"
    "S92_CD3_D12_SFEM__1312024_S92_ME_L001"
    "S93_CP3_D7_SFEM__1312024_S93_ME_L001"
    "S94_CP3_D14_SFEM__1312024_S94_ME_L001"
)

# Loop through each sample and align with BWA MEM
for sample in "${samples[@]}"; do
    TRIMMED_R1="${input_dir}/${sample}_R1_trimmed.fastq.gz"
    TRIMMED_R2="${input_dir}/${sample}_R2_trimmed.fastq.gz"
    BAM_OUTPUT="${output_dir}/${sample}.bam"

    # Align with BWA MEM and convert SAM to BAM
    if [[ -e $TRIMMED_R1 && -e $TRIMMED_R2 ]]; then
        bwa mem -M -t ${SLURM_CPUS_PER_TASK} $HUMAN_GENOME $TRIMMED_R1 $TRIMMED_R2 | samtools view -bS - > $BAM_OUTPUT
        echo "Alignment completed for ${sample}"
    else
        echo "Trimmed files for $sample are missing!"
    fi
done

echo "All samples processed."
