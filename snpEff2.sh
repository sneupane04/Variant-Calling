#!/bin/bash
#SBATCH --job-name=snpEff_annotation
#SBATCH --mail-user=surendra.neupane@moffitt.org
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --time=24:00:00

# Load the snpEff module
module load snpEff/5.2a

# Define the directory containing VCF files
VCF_DIR="/share/lab_padron/CICPT-4562-Padron_Traci-QIAseq-Dec-01-2023_Data/fastq/fastqtraci/bothsam"


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

# Copy the snpEff configuration file to the current directory
cp $EBROOTSNPEFF/snpEff.config ./

# Iterate over each sample
for sample in "${samples[@]}"; do
    INPUT_VCF="${VCF_DIR}/${sample}_variants_bcftools.vcf"
    ANNOTATED_VCF="${VCF_DIR}/${sample}_annotated_variants_bcftools.vcf"

    if [[ -e $INPUT_VCF ]]; then
        # Annotate the VCF file with snpEff
        echo "Annotating VCF file with snpEff for $sample"
        java -Xmx4g -jar $EBROOTSNPEFF/snpEff.jar -config snpEff.config GRCh38.p14 "$INPUT_VCF" > "$ANNOTATED_VCF"
    else
        echo "VCF file for $sample is missing!"
    fi
done
