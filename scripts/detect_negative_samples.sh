#!/bin/bash

# detect_negative_samples.sh - Identify samples with insufficient coverage after first pass alignment

# Coverage threshold (average coverage below this is considered negative)
COVERAGE_THRESHOLD=${1:-50}  # Default to 50x coverage

echo "Checking samples for negative status (threshold: ${COVERAGE_THRESHOLD}x average coverage)"

# Working in first_pass results
FIRST_PASS_DIR="results/first_pass"
NEGATIVE_DIR="results/negative_samples"
TRIMMED_DIR="results/trimmed"  # New centralized location for trimmed reads

# Track identified negative samples
NEGATIVE_SAMPLES_FILE="$NEGATIVE_DIR/negative_samples.txt"
mkdir -p "$NEGATIVE_DIR"
> "$NEGATIVE_SAMPLES_FILE"  # Create empty file

echo "DEBUG: Negative samples will be listed in $NEGATIVE_SAMPLES_FILE"

# Process each sample in first_pass directory
for sample_dir in "$FIRST_PASS_DIR"/*; do
    if [ ! -d "$sample_dir" ]; then
        continue
    fi
    
    sample=$(basename "$sample_dir")
    
    # Skip if not a sample directory
    if [[ "$sample" == "plots" ]]; then
        continue
    fi
    
    echo "Analyzing coverage for sample: $sample"
    
    # Track if this sample is negative across all segments
    is_negative=true
    sample_has_bam=false
    avg_coverage_info=""
    
    # Check each segment
    for segment in L_segment M_segment S_segment; do
        bam_file="$FIRST_PASS_DIR/$sample/alignment/$segment/${sample}.bam"
        
        if [ -f "$bam_file" ]; then
            sample_has_bam=true
            
            # Calculate average coverage using samtools
            avg_coverage=$(samtools depth -a "$bam_file" | awk '{sum+=$3} END {if(NR>0) print sum/NR; else print 0}')
            avg_coverage_rounded=$(printf "%.1f" "$avg_coverage")
            
            echo "  $segment: Average coverage = ${avg_coverage_rounded}x"
            avg_coverage_info="${avg_coverage_info}  $segment: ${avg_coverage_rounded}x\n"
            
            # If any segment has coverage above threshold, sample is not negative
            if (( $(echo "$avg_coverage > $COVERAGE_THRESHOLD" | bc -l) )); then
                echo "  $segment has coverage above threshold - marking as POSITIVE"
                is_negative=false
            fi
        fi
    done
    
    # Only process if we have alignment files
    if [ "$sample_has_bam" = true ]; then
        if [ "$is_negative" = true ]; then
            echo "  RESULT: NEGATIVE sample detected ($sample)"
            echo "$sample" >> "$NEGATIVE_SAMPLES_FILE"
            echo -e "$avg_coverage_info"
            
            # Create directory for negative sample
            mkdir -p "$NEGATIVE_DIR/$sample"
            
            # Copy alignment files from first pass
            echo "  Copying alignment files to negative samples directory"
            cp -r "$FIRST_PASS_DIR/$sample/alignment" "$NEGATIVE_DIR/$sample/"
            
            # Also copy consensus files if they exist
            if [ -d "$FIRST_PASS_DIR/$sample/consensus" ]; then
                mkdir -p "$NEGATIVE_DIR/$sample/consensus"
                cp -r "$FIRST_PASS_DIR/$sample/consensus" "$NEGATIVE_DIR/$sample/"
            fi
            
            # Also copy variant calls if they exist
            if [ -d "$FIRST_PASS_DIR/$sample/variants" ]; then
                mkdir -p "$NEGATIVE_DIR/$sample/variants"
                cp -r "$FIRST_PASS_DIR/$sample/variants" "$NEGATIVE_DIR/$sample/"
            fi
            
            # Add a marker file to ensure this sample is treated as negative
            touch "$NEGATIVE_DIR/$sample/.is_negative"
        else
            echo "  RESULT: POSITIVE sample ($sample) - will proceed to second pass"
        fi
    fi
done

# Count negative samples
neg_count=$(wc -l < "$NEGATIVE_SAMPLES_FILE")
echo "Detection complete. Found $neg_count negative samples."

# Print the list of negative samples for debugging
if [ "$neg_count" -gt 0 ]; then
    echo "List of negative samples:"
    cat "$NEGATIVE_SAMPLES_FILE"
fi 