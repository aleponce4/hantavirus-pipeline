#!/bin/bash
# Filter variants for primer mismatch detection
# This script is used by the main variant calling script

sample=$1
segment=$2
results_dir=$3
min_quality=${4:-20}    # Default to 20 if not provided
min_freq=${5:-0.03}     # Default to 0.03 if not provided

if [ -z "$sample" ] || [ -z "$segment" ] || [ -z "$results_dir" ]; then
    echo "Usage: $0 <sample> <segment> <results_dir> [min_quality] [min_freq]"
    exit 1
fi

# Define paths
variant_tsv="$results_dir/$sample/variants/$segment/${sample}.tsv"
filtered_tsv="$results_dir/$sample/variants/$segment/${sample}.filtered.tsv"

# Check if input exists
if [ ! -f "$variant_tsv" ]; then
    echo "Error: Variant file not found: $variant_tsv"
    exit 1
fi

# Filter variants based on quality and frequency
echo "  Filtering variants with quality >= $min_quality and frequency >= $min_freq"
cat "$variant_tsv" | awk -v qual="$min_quality" -v freq="$min_freq" 'NR==1 || ($10 >= qual && $11 >= freq)' > "$filtered_tsv"

echo "  Filtered variants saved to: $filtered_tsv"
echo "  $(grep -v "REGION" "$filtered_tsv" | wc -l) variants passed filtering." 