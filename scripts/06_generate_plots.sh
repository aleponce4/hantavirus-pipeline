#!/bin/bash

# Script to generate alignment and coverage visualizations

echo "Generating visualizations"

# Create plots directory if it doesn't exist
mkdir -p results/plots

# Check for MAFFT and install if needed
if ! command -v mafft &> /dev/null; then
    echo "MAFFT not found. Installing MAFFT..."
    conda install -y -c bioconda mafft
fi

# 1. Create alignment visualizations
echo "Creating alignment visualizations..."
python3 scripts/plotting/plot_metaconsensus_comparison.py

# 2. Create coverage plots for all samples
echo "Creating coverage plots..."
python3 scripts/plotting/plot_coverage.py

echo "Visualizations complete"
echo "Results are available in results/plots/" 