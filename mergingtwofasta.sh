#!/bin/bash
#SBATCH --job-name=STAR
#SBATCH --mail-user=surendra.neupane@moffitt.org
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

date;hostname;pwd

### MERGE TWO GENOMES (FASTA FORMAT)

# specify genome path
GENOME_A="/share/lab_padron/CICPT-4562-Padron_Traci-QIAseq-Dec-01-2023_Data/fastq/fastqtraci/GRCh38.p14.genome.fa"
GENOME_B="/share/lab_padron/CICPT-4562-Padron_Traci-QIAseq-Dec-01-2023_Data/fastq/fastqtraci/GRCm39.genome.fa"
OUTPUT_GENOME="/share/lab_padron/CICPT-4562-Padron_Traci-QIAseq-Dec-01-2023_Data/fastq/fastqtraci/merged_ICRG.fasta"

# generate a copy of genome_b with modified chromosome names
cp "$GENOME_B" "${GENOME_B}_edited.fasta"
if sed -i 's/^>chr/>m.chr/g' "${GENOME_B}_edited.fasta"; then
    echo "Chromosome names in $GENOME_B have been successfully modified."
else
    echo "Error modifying chromosome names in $GENOME_B."
    exit 1
fi

# combine the two fasta files
if cat "$GENOME_A" "${GENOME_B}_edited.fasta" > "$OUTPUT_GENOME"; then
    echo "Genomes have been successfully merged into $OUTPUT_GENOME."
else
    echo "Error merging genomes."
    exit 1
fi

