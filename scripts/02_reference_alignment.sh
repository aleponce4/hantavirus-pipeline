#!/bin/bash

# Reference alignment script

sample=$1
bwa_threads=$2
samtools_threads=$3

# Check if this is the first pass
if [ "$FIRST_PASS" = "true" ]; then
    echo "Aligning sample: $sample with BWA threads: $bwa_threads, Samtools threads: $samtools_threads (FIRST PASS - using original reference)"
    RESULTS_DIR="results/first_pass"
else
    echo "Aligning sample: $sample with BWA threads: $bwa_threads, Samtools threads: $samtools_threads (SECOND PASS - using metaconsensus)"
    RESULTS_DIR="results/second_pass"
fi

# Input files - use environment variable if provided or default location
TRIMMED_DIR=${TRIMMED_DIR:-"results/$sample/trimmed"}
R1_PAIRED="$TRIMMED_DIR/${sample}_R1_paired.fastq.gz"
R2_PAIRED="$TRIMMED_DIR/${sample}_R2_paired.fastq.gz"

# Check if files exist
if [ ! -f "$R1_PAIRED" ] || [ ! -f "$R2_PAIRED" ]; then
    echo "ERROR: Trimmed files not found at $TRIMMED_DIR"
    echo "Looked for: $R1_PAIRED and $R2_PAIRED"
    exit 1
fi

# Output directory
OUT_DIR="$RESULTS_DIR/$sample/alignment"
mkdir -p $OUT_DIR

# Check for segments
for segment in L_segment M_segment S_segment; do
    segment_dir="data/references/$segment"
    
    # Create segment output directory
    seg_out_dir="$OUT_DIR/$segment"
    mkdir -p $seg_out_dir
    
    # Check if there are any files in the segment directory
    if [ -d "$segment_dir" ] && [ "$(ls -A $segment_dir)" ]; then
        # Reference selection strategy
        if [ "$FIRST_PASS" = "true" ]; then
            # For first pass, always use the original reference
            reference=$(ls $segment_dir/*.fasta | grep -v "metaconsensus.fasta" | head -n 1)
            echo "FIRST PASS: Using original reference for $segment: $reference"
        else
            # For second pass, always use the metaconsensus reference
            metaconsensus="data/references/$segment/metaconsensus.fasta"
            if [ -f "$metaconsensus" ]; then
                reference=$metaconsensus
                echo "SECOND PASS: Using metaconsensus reference for $segment: $reference"
            else
                echo "ERROR: Metaconsensus reference not found for $segment. Cannot proceed with second pass."
                exit 1
            fi
        fi
        
        if [ -n "$reference" ]; then
            # Index reference
            bwa index $reference
            samtools faidx $reference
            
            # Align reads to reference using multiple threads
            # BWA for alignment, samtools for sorting with separate thread counts
            bwa mem -t $bwa_threads -k 15 -B 3 $reference $R1_PAIRED $R2_PAIRED | \
                samtools sort -@ $samtools_threads -o "$seg_out_dir/${sample}.bam"
            
            # Index BAM file
            samtools index "$seg_out_dir/${sample}.bam"
            
            echo "Alignment completed for sample $sample, segment $segment"
        fi
    fi
done

echo "Alignment completed for sample: $sample" 