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

echo "Pipeline finished. Results are in the 'results' directory."

# Run primer evaluation
echo "Running primer evaluation..."
bash ./scripts/07_evaluate_primers.sh
echo "Primer evaluation finished. Results are in the 'results/primer_evaluation' directory." 