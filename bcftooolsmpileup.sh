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
BAM_DIR="/share/lab_padron/CICPT-4562-Padron_Traci-QIAseq-Dec-01-2023_Data/fastq/fastqtraci/bothsam"

# Define the samples array
samples=(
    "CICPT_4562_2666_S26-11242023_S11_ME_L001"
    "CICPT_4562_2664_S25-11242023_S4_ME_L001"
    "CICPT_4562_2662_S24-11242023_S52_ME_L001"
    "CICPT_4562_2660_S23-11242023_S45_ME_L001"
    "CICPT_4562_2658_S22-11242023_S38_ME_L001"
    "CICPT_4562_2656_S21-11242023_S31_ME_L001"
    "CICPT_4562_2654_S20-11242023_S24_ME_L001"
    "CICPT_4562_2652_S19-11242023_S17_ME_L001"
    "CICPT_4562_2630_S27-11242023_S18_ME_L001"
    "CICPT_4562_2614_S18-11242023_S10_ME_L001"
    "CICPT_4562_2612_S17-11242023_S3_ME_L001"
    "CICPT_4562_2610_S16-11242023_S51_ME_L001"
    "CICPT_4562_2608_S15-11242023_S44_ME_L001"
    "CICPT_4562_2606_S14-11242023_S37_ME_L001"
    "CICPT_4562_2604_S13-11242023_S30_ME_L001"
    "CICPT_4562_2602_S12-11242023_S23_ME_L001"
    "CICPT_4562_2600_S11-11242023_S16_ME_L001"
    "CICPT_4562_2572_S31-11242023_S46_ME_L001"
    "CICPT_4562_2562_S10-11242023_S9_ME_L001"
    "CICPT_4562_2558_S9-11242023_S2_ME_L001"
    "CICPT_4562_2534_S30-11242023_S39_ME_L001"
    "CICPT_4562_2526_S29-11242023_S32_ME_L001"
    "CICPT_4562_2498_S32-11242023_S53_ME_L001"
    "CICPT_4562_2488_S28-11242023_S25_ME_L001"
    "CICPT_4562_2334_S8-11242023_S50_ME_L001"
    "CICPT_4562_2332_S7-11242023_S43_ME_L001"
    "CICPT_4562_2330_S6-11242023_S36_ME_L001"
    "CICPT_4562_2328_S5-11242023_S29_ME_L001"
    "CICPT_4562_2326_S4-11242023_S22_ME_L001"
    "CICPT_4562_2324_S3-11242023_S15_ME_L001"
    "CICPT_4562_2322_S2-11242023_S8_ME_L001"
    "CICPT_4562_2320_S1-11242023_S1_ME_L001"
)

# Iterate over each sample
for sample in "${samples[@]}"; do
    INPUT_BAM="${BAM_DIR}/${sample}_genome_human_only_dedup.bam"
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

