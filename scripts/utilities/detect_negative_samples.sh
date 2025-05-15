#!/bin/bash
# Script to detect negative samples based on low coverage

# Source the configuration file
source ./config.sh

# Coverage threshold (from config or from command line)
COVERAGE_THRESHOLD=${1:-$NEGATIVE_SAMPLE_THRESHOLD}

echo "Checking samples for negative status (threshold: ${COVERAGE_THRESHOLD}x average coverage)"

# Output directory
mkdir -p results/negative_samples

# Create a file to store negative sample IDs
> results/negative_samples/negative_samples.txt

# Create a directory for negative sample coverage plots
mkdir -p results/plots/negative_samples

# Check all samples in results/first_pass (using first pass results for classification)
for sample_dir in results/first_pass/*; do
    if [ ! -d "$sample_dir" ]; then
        continue
    fi
    
    sample=$(basename "$sample_dir")
    echo "Checking coverage for sample: $sample"
    
    # Flag to track if sample is negative
    is_negative=true
    
    # Check each segment
    for segment in S_segment M_segment L_segment; do
        bam_file="results/first_pass/$sample/alignment/$segment/${sample}.bam"
        
        if [ ! -f "$bam_file" ]; then
            echo "  No BAM file found for $segment, skipping"
            continue
        fi
        
        # Generate coverage data
        coverage_file="results/negative_samples/${sample}_${segment}_coverage.txt"
        samtools depth -a "$bam_file" > "$coverage_file"
        
        # Check if coverage file is empty
        if [ ! -s "$coverage_file" ]; then
            echo "  No coverage data for $segment - marking as NEGATIVE"
            continue
        fi
        
        # Calculate average coverage
        avg_coverage=$(awk '{sum+=$3} END {print sum/NR}' "$coverage_file")
        printf "  %s average coverage: %.2f\n" "$segment" "$avg_coverage"
        
        # If any segment has coverage above threshold, sample is not negative
        if (( $(echo "$avg_coverage > $COVERAGE_THRESHOLD" | bc -l) )); then
            echo "  $segment has coverage above threshold - marking as POSITIVE"
            is_negative=false
        fi
        
        # Generate coverage plot for potential negative samples
        if (( $(echo "$avg_coverage < $COVERAGE_THRESHOLD * 2" | bc -l) )); then
            echo "  Generating coverage plot for low-coverage segment"
            plot_file="results/plots/negative_samples/${sample}_${segment}_coverage.png"
            python scripts/plotting/plot_coverage.py --bam "$bam_file" --output "$plot_file" --title "$sample - $segment Coverage" --low_coverage $COVERAGE_THRESHOLD
        fi
    done
    
    # If all segments have low coverage, mark as negative
    if [ "$is_negative" = true ]; then
        echo "$sample" >> results/negative_samples/negative_samples.txt
        echo "Sample $sample classified as NEGATIVE"
        
        # Copy any variant calls and consensus sequences to negative_samples directory
        for segment in S_segment M_segment L_segment; do
            # Copy variant files if they exist
            variant_dir="results/first_pass/$sample/variants/$segment"
            if [ -d "$variant_dir" ]; then
                mkdir -p "results/negative_samples/$sample/variants/$segment"
                cp -f "$variant_dir"/*.tsv "results/negative_samples/$sample/variants/$segment/" 2>/dev/null || true
                cp -f "$variant_dir"/*.vcf.gz "results/negative_samples/$sample/variants/$segment/" 2>/dev/null || true
                cp -f "$variant_dir"/*.vcf.gz.tbi "results/negative_samples/$sample/variants/$segment/" 2>/dev/null || true
            fi
            
            # Copy consensus files if they exist
            consensus_dir="results/first_pass/$sample/consensus/$segment"
            if [ -d "$consensus_dir" ]; then
                mkdir -p "results/negative_samples/$sample/consensus/$segment"
                cp -f "$consensus_dir"/*.fasta "results/negative_samples/$sample/consensus/$segment/" 2>/dev/null || true
            fi
        done
    else
        echo "Sample $sample classified as POSITIVE"
    fi
done

# Count and report negative samples
negative_count=$(wc -l < results/negative_samples/negative_samples.txt)
echo "Identified $negative_count negative samples (average coverage < ${COVERAGE_THRESHOLD}x)"
if [ "$negative_count" -gt 0 ]; then
    echo "Negative sample IDs:"
    cat results/negative_samples/negative_samples.txt
fi 