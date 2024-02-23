#!/bin/bash
#SBATCH --job-name=haplotype_caller
#SBATCH --mail-user=surendra.neupane@moffitt.org
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --time=24:00:00

# Load modules
module load GATK/4.2.0.0-GCCcore-10.2.0-Java-11
module load picard/2.25.5-Java-13
module load SAMtools/1.16.1-GCC-11.3.0
module load BCFtools/1.12-GCC-10.2.0

# Define the reference genome path
REF_GENOME="/share/lab_padron/CICPT-4562-Padron_Traci-QIAseq-Dec-01-2023_Data/fastq/fastqtraci/GRCh38.p14.genome.fa"

# Check if FASTA index and dictionary files exist, if not create them
if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "Creating FASTA index for ${REF_GENOME}"
    samtools faidx "$REF_GENOME"
fi
if [ ! -f "${REF_GENOME%.fa}.dict" ]; then
    echo "Creating FASTA dictionary for ${REF_GENOME}"
    java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R="$REF_GENOME" O="${REF_GENOME%.fa}.dict"
fi

# Define the directory containing BAM files
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

# Iterate over each sample
for sample in "${samples[@]}"; do
    INPUT_BAM="${BAM_DIR}/${sample}_dedup.bam"
    OUTPUT_VCF="${BAM_DIR}/${sample}_variants_bcftools.vcf"

    if [[ -e $INPUT_BAM ]]; then
        # Index the BAM file
        echo "Indexing $INPUT_BAM"
        samtools index "$INPUT_BAM"

        # Call variants with SAMtools and BcfTools
        echo "Running SAMtools mpileup and BcfTools view on $INPUT_BAM"
        bcftools mpileup -f "$REF_GENOME" "$INPUT_BAM" | bcftools call -mv -Ov > "$OUTPUT_VCF"
    else
        echo "Deduplicated BAM file for $sample is missing!"
    fi
done

