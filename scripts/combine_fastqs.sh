#!/bin/bash
#SBATCH --job-name=combine_fastqs
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# Hardcoded variables
INPUT_R1=("sample1_R1.fq.gz" "sample2_R1.fq.gz" "sample3_R1.fq.gz")
OUTPUT_BASE="combined_samplename"  # Base name for the output files

# Generate output file names dynamically
OUTPUT_R1="${OUTPUT_BASE}_R1.fastq.gz"
OUTPUT_R2="${OUTPUT_BASE}_R2.fastq.gz"

# Generate the list of R2 file names from R1
INPUT_R2=()
for fq in "${INPUT_R1[@]}"; do
  INPUT_R2+=("${fq/_R1/_R2}")
done

# Combining R1 files
echo "Combining FASTQ files for R1..."
zcat "${INPUT_R1[@]}" | gzip > "$OUTPUT_R1"

# Combining R2 files
echo "Combining FASTQ files for R2..."
zcat "${INPUT_R2[@]}" | gzip > "$OUTPUT_R2"

echo "FASTQ files combined successfully!"
