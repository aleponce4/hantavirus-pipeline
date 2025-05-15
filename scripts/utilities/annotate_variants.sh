#!/bin/bash
# Script to normalize and annotate VCF files with SnpEff

# Parameters
sample=$1
segment=$2
results_dir=$3

# Source configuration file
source ./config.sh

if [ -z "$sample" ] || [ -z "$segment" ] || [ -z "$results_dir" ]; then
    echo "Usage: $0 <sample> <segment> <results_dir>"
    exit 1
fi

# Check if SnpEff is installed
if [ ! -f "$SNPEFF_DIR/snpEff.jar" ]; then
    echo "SnpEff not found. Installing..."
    bash scripts/utilities/setup_snpeff.sh
fi

# Define file paths
segment_short=$(echo "$segment" | cut -d'_' -f1)
vcf_file="$results_dir/$sample/variants/$segment/${sample}.vcf.gz"
unzipped_vcf="$results_dir/$sample/variants/$segment/${sample}.vcf"
annotated_vcf="$results_dir/$sample/variants/$segment/${sample}.annotated.vcf"

# Check if the VCF file exists
if [ ! -f "$vcf_file" ]; then
    echo "Error: VCF file not found: $vcf_file"
    exit 1
fi

# Unzip the VCF file if it's compressed
if [ -f "$vcf_file" ] && [ ! -f "$unzipped_vcf" ]; then
    echo "Unzipping VCF file..."
    gunzip -c "$vcf_file" > "$unzipped_vcf"
fi

# Prepare for annotation
echo "Normalizing and annotating variants for $sample - $segment"

# Change to SnpEff directory
cd "$SNPEFF_DIR"

# Step 1: Normalize and annotate the VCF file with SnpEff in a single step
# Use correct capitalization for SnpEff flags
java -Xmx${SNPEFF_MEMORY} -jar snpEff.jar -noDownload -dataDir ./data -v \
    "hanta_${segment_short}_segment" "../../$unzipped_vcf" > "../../$annotated_vcf"

# Go back to original directory
cd - > /dev/null

# Step 2: Compress and index the annotated VCF
bgzip -f -c "$annotated_vcf" > "$annotated_vcf.gz"
tabix -p vcf "$annotated_vcf.gz"

# Step 3: Replace the original VCF with the annotated one
mv "$annotated_vcf.gz" "$vcf_file"
mv "$annotated_vcf.gz.tbi" "$vcf_file.tbi"

# Clean up
rm -f "$unzipped_vcf" "$annotated_vcf"

echo "Variant annotation completed for $sample - $segment" 