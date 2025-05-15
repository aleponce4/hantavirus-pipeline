#!/bin/bash

# Preprocessing script - Trims raw reads using Trim Galore

sample=$1
threads=$2

# Source the configuration file
source ./config.sh

echo "Preprocessing sample: $sample"
echo "Using $threads threads"

# Input files
R1="data/raw_reads/${sample}_R1.fastq.gz"
R2="data/raw_reads/${sample}_R2.fastq.gz"

# Output files - use environment variable if provided or default location
OUT_DIR=${OUT_DIR:-"results/$sample/trimmed"}
mkdir -p $OUT_DIR

# Print conda environment info
echo "Using conda environment: $CONDA_DEFAULT_ENV"
echo "CONDA_PREFIX: $CONDA_PREFIX"
echo "PATH: $PATH"

# Verify trim_galore is available
if ! command -v trim_galore &> /dev/null; then
    echo "Error: trim_galore command not found"
    exit 1
fi

# Run Trim Galore with multithreading
echo "Running Trim Galore with $threads threads..."
trim_galore --paired --quality $TRIM_QUALITY --stringency $ADAPTER_STRINGENCY --length $MIN_READ_LENGTH \
    --output_dir $OUT_DIR \
    --basename ${sample} \
    --cores $threads \
    $R1 $R2

# Rename output files to match expected names in later steps
mv "$OUT_DIR/${sample}_val_1.fq.gz" "$OUT_DIR/${sample}_R1_paired.fastq.gz"
mv "$OUT_DIR/${sample}_val_2.fq.gz" "$OUT_DIR/${sample}_R2_paired.fastq.gz"

echo "Quality trimming completed for sample: $sample" 