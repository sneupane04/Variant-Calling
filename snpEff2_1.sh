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
VCF_DIR=""/share/lab_padron/CICPT_4628_AgilentPanel-February-2024/merged""


# Define the samples array
samples=("S96_Enzo-Positive-Control--20ng_Enzo-Positive-Control--20ng_1312024_S6S96_ME"
    "S87_CP1_D7_SFEM__1312024_S87S3_ME"
    "S90_CD2_D7_SFEM__1312024_S4S90_ME"
    "S88_CP1_D14_SFEM__1312024_S88_ME_L001"
    "S89_CP2_D7_SFEM__1312024_S89_ME_L001"
    "S91_CD3_D7_SFEM__1312024_S91_ME_L001"
    "S92_CD3_D12_SFEM__1312024_S92_ME_L001"
    "S93_CP3_D7_SFEM__1312024_S93_ME_L001"
    "S94_CP3_D14_SFEM__1312024_S94_ME_L001"  
)

# Copy the snpEff configuration file to the current directory
cp $EBROOTSNPEFF/snpEff.config ./

# Iterate over each sample
for sample in "${samples[@]}"; do
    INPUT_VCF="${VCF_DIR}/${sample}_variants_GATK.vcf"
    ANNOTATED_VCF="${VCF_DIR}/${sample}_annotated_variants_GATK.vcf"

    if [[ -e $INPUT_VCF ]]; then
        # Annotate the VCF file with snpEff
        echo "Annotating VCF file with snpEff for $sample"
        java -Xmx4g -jar $EBROOTSNPEFF/snpEff.jar -config snpEff.config GRCh38.p14 "$INPUT_VCF" > "$ANNOTATED_VCF"
    else
        echo "VCF file for $sample is missing!"
    fi
done
