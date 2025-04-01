#!/bin/bash

# Script to generate alignment and coverage visualizations

echo "Generating visualizations"

# Create plots directory if it doesn't exist
mkdir -p results/plots
mkdir -p results/plots/negative_samples

# Check for MAFFT and install if needed
if ! command -v mafft &> /dev/null; then
    echo "MAFFT not found. Installing MAFFT..."
    conda install -y -c bioconda mafft
fi

# 1. Create alignment visualizations for positive samples
echo "Creating alignment visualizations..."
python3 scripts/plotting/plot_metaconsensus_comparison.py

# 2. Create coverage plots for positive samples
echo "Creating coverage plots for positive samples..."
python3 scripts/plotting/plot_coverage.py --results_dir "results/second_pass"

# 3. Create coverage plots for negative samples if any exist
if [ -d "results/negative_samples" ] && [ "$(ls -A results/negative_samples 2>/dev/null)" ]; then
    echo "Creating coverage plots for negative samples..."
    for segment in L_segment M_segment S_segment; do
        python3 scripts/plotting/plot_coverage.py \
            --results_dir "results/negative_samples" \
            --segment "$segment" \
            --output "results/plots/negative_samples/${segment}_coverage.png" \
            --max_coverage 100  # Lower max coverage to see details better
    done
fi

echo "Visualizations complete"
echo "Results are available in results/plots/" 