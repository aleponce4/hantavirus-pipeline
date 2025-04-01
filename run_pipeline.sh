#!/bin/bash

# Simple wrapper script to run the viral sequencing pipeline

# Create logs directory
mkdir -p logs

echo "Starting viral sequencing data processing pipeline ($(date))"

# Ensure proper conda activation
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate De_Novo_pipeline

# Verify conda environment is active
if [[ "${CONDA_DEFAULT_ENV}" != "De_Novo_pipeline" ]]; then
    echo "Error: Failed to activate De_Novo_pipeline conda environment"
    exit 1
fi
echo "Successfully activated conda environment: ${CONDA_DEFAULT_ENV}"

# Run the main pipeline
bash ./pipeline.sh

# Check for negative samples
if [ -f "results/negative_samples/negative_samples.txt" ]; then
    neg_count=$(wc -l < "results/negative_samples/negative_samples.txt")
    if [ "$neg_count" -gt 0 ]; then
        echo "Detected $neg_count negative samples (see results/negative_samples/)"
        echo "Coverage plots for negative samples are in results/plots/negative_samples/"
    fi
fi

# Clean up any sample directories in the root of results that are not in the proper structure
echo "Cleaning up results directory structure..."
for item in results/*; do
    basename=$(basename "$item")
    # Only remove directories that appear to be sample directories
    if [ -d "$item" ] && [ "$basename" != "first_pass" ] && [ "$basename" != "second_pass" ] && \
       [ "$basename" != "negative_samples" ] && [ "$basename" != "plots" ] && \
       [ "$basename" != "primer_evaluation" ] && [ "$basename" != "trimmed" ]; then
        echo "  Removing leftover directory: $item (data preserved in organized structure)"
        rm -rf "$item"
    fi
done

echo "Pipeline finished. Results are organized in the 'results' directory with the following structure:"
echo "  - results/first_pass/: First pass processing results"
echo "  - results/second_pass/: Second pass processing results (positive samples only)"
echo "  - results/negative_samples/: Negative sample data"
echo "  - results/plots/: Coverage and alignment plots"
echo "  - results/trimmed/: Trimmed reads (used by both passes)"

# Run primer evaluation
echo "Running primer evaluation..."
bash ./scripts/07_evaluate_primers.sh
echo "Primer evaluation finished. Results are in the 'results/primer_evaluation' directory." 