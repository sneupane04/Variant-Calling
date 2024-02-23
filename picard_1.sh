#!/bin/bash
#SBATCH --job-name=remove_duplicates
#SBATCH --mail-user=surendra.neupane@moffitt.org
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --time=24:00:00

# Load necessary modules (adjust according to your environment)
# module load picard or set PICARD_JAR to the path of the Picard jar file

module load picard/2.25.5-Java-13
module load SAMtools/1.16.1-GCC-11.3.0

# Define the directory containing the BAM files
BAM_DIR="/share/lab_padron/CICPT_4628_AgilentPanel-February-2024/merged"

# Define the samples array
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

# Iterate over each sample and remove duplicates
for sample in "${samples[@]}"; do
    INPUT_BAM="${BAM_DIR}/${sample}.bam"
    SORTED_BAM="${BAM_DIR}/${sample}_sorted.bam"
    OUTPUT_BAM="${BAM_DIR}/${sample}_dedup.bam"
    METRICS_FILE="${BAM_DIR}/${sample}_dedup_metrics.txt"

    if [[ -e $INPUT_BAM ]]; then
        # Sort the BAM file by coordinates
        samtools sort -o $SORTED_BAM $INPUT_BAM

        # Remove duplicates using Picard
        java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
            I=$SORTED_BAM \
            O=$OUTPUT_BAM \
            M=$METRICS_FILE \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=SILENT
    else
        echo "BAM file for $sample is missing!"
    fi
done
