#!/bin/bash

# Preprocessing script - Trims raw reads using Trim Galore and removes primers using cutadapt

sample=$1
threads=$2

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
trim_galore --paired --quality 20 --stringency 3 --length 36 \
    --output_dir $OUT_DIR \
    --basename ${sample} \
    --cores $threads \
    $R1 $R2

# Rename output files to match expected names in later steps
mv "$OUT_DIR/${sample}_val_1.fq.gz" "$OUT_DIR/${sample}_R1_paired.fastq.gz"
mv "$OUT_DIR/${sample}_val_2.fq.gz" "$OUT_DIR/${sample}_R2_paired.fastq.gz"

echo "Quality trimming completed for sample: $sample"

# Verify cutadapt is available
if ! command -v cutadapt &> /dev/null; then
    echo "Error: cutadapt command not found"
    exit 1
fi

# Define temporary files for the primer trimming step
TEMP_R1="$OUT_DIR/${sample}_R1_temp.fastq.gz"
TEMP_R2="$OUT_DIR/${sample}_R2_temp.fastq.gz"
TRIMMED_R1="$OUT_DIR/${sample}_R1_paired.fastq.gz"
TRIMMED_R2="$OUT_DIR/${sample}_R2_paired.fastq.gz"

# Define primer lists for each segment
# S segment primers (forward and reverse)
S_PRIMERS=()
# M segment primers (forward and reverse)
M_PRIMERS=()
# L segment primers (forward and reverse)
L_PRIMERS=()

# Extract primers from the primers file
echo "Extracting segment-specific primers..."
while IFS=, read -r name region sequence rest; do
    # Skip header line
    if [[ $name == "F Name" ]]; then
        continue
    fi
    
    # Clean up sequence (remove spaces)
    sequence=$(echo "$sequence" | tr -d ' ')
    
    # Add to appropriate segment list based on primer name
    if [[ $name == SF* || $name == SR* ]]; then
        S_PRIMERS+=("$sequence")
    elif [[ $name == MF* || $name == MR* ]]; then
        M_PRIMERS+=("$sequence")
    elif [[ $name == LF* || $name == LR* ]]; then
        L_PRIMERS+=("$sequence")
    fi
done < "data/primers/Primers.csv"

# Process each segment
for segment in "S_segment" "M_segment" "L_segment"; do
    echo "Processing $segment primers..."
    
    # Set the primers list based on segment
    if [[ $segment == "S_segment" ]]; then
        PRIMERS=("${S_PRIMERS[@]}")
    elif [[ $segment == "M_segment" ]]; then
        PRIMERS=("${M_PRIMERS[@]}")
    elif [[ $segment == "L_segment" ]]; then
        PRIMERS=("${L_PRIMERS[@]}")
    fi
    
    # Skip if no primers for this segment
    if [ ${#PRIMERS[@]} -eq 0 ]; then
        echo "No primers found for $segment, skipping."
        continue
    fi
    
    # Build cutadapt command with all primers for this segment
    CUTADAPT_ARGS=""
    for primer in "${PRIMERS[@]}"; do
        # Add primer to be trimmed from both ends of both reads
        CUTADAPT_ARGS+=" -a $primer -A $primer -g $primer -G $primer"
    done
    
    # Run cutadapt for this segment
    echo "Trimming $segment primers from reads..."
    cutadapt $CUTADAPT_ARGS \
        -j $threads \
        -o $TEMP_R1 -p $TEMP_R2 \
        $TRIMMED_R1 $TRIMMED_R2 \
        > "$OUT_DIR/${sample}_${segment}_primer_trimming.log"
    
    # Extract summary statistics from log - using the standard cutadapt summary format
    removed=$(grep "Reads with adapters:" "$OUT_DIR/${sample}_${segment}_primer_trimming.log" | head -1 | sed 's/Reads with adapters: //')
    
    # If empty, use a default value
    if [ -z "$removed" ]; then
        removed="0 (0.0%)"
    fi
    
    echo "  $segment primer trimming complete: $removed reads had primers removed"
    
    # Update the trimmed files for the next segment
    mv $TEMP_R1 $TRIMMED_R1
    mv $TEMP_R2 $TRIMMED_R2
done

echo "Preprocessing completed for sample: $sample" 