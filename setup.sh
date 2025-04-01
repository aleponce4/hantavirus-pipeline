#!/bin/bash

# Hantavirus Pipeline Setup Script
# This script sets up the conda environment and verifies required dependencies

set -e  # Exit immediately if a command exits with a non-zero status

echo "Setting up Hantavirus Pipeline Environment..."

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "Error: conda is not installed or not in PATH"
    echo "Please install Miniconda or Anaconda first: https://docs.conda.io/projects/conda/en/latest/user-guide/install/"
    exit 1
fi

# Check if environment.yml exists
if [ ! -f environment.yml ]; then
    echo "Error: environment.yml file not found"
    exit 1
fi

# Create conda environment from YAML file
echo "Creating conda environment 'De_Novo_pipeline' from environment.yml..."
conda env create -f environment.yml

echo "Verifying installation of key dependencies..."

# Activate environment and check for key dependencies
eval "$(conda shell.bash hook)"
conda activate De_Novo_pipeline

# Check key dependencies
DEPENDENCIES=("trim_galore" "samtools" "bwa" "cutadapt" "python" "seqkit" "lofreq")
MISSING=()

for dep in "${DEPENDENCIES[@]}"; do
    if ! command -v $dep &> /dev/null; then
        MISSING+=($dep)
    fi
done

if [ ${#MISSING[@]} -eq 0 ]; then
    echo "All dependencies are properly installed"
else
    echo "Warning: The following dependencies could not be found in the environment:"
    for missing in "${MISSING[@]}"; do
        echo "  - $missing"
    done
    echo "You may need to install them manually or fix your environment"
fi

# Create necessary directories if they don't exist
echo "Creating required directories..."
mkdir -p data/raw_reads
mkdir -p data/references
mkdir -p data/primers
mkdir -p results

echo ""
echo "Setup complete!"
echo ""
echo "To use the pipeline:"
echo "1. Activate the environment:   conda activate De_Novo_pipeline"
echo "2. Place raw reads in:         data/raw_reads/"
echo "3. Place references in:        data/references/"
echo "4. Run the pipeline:           bash run_pipeline.sh [sample_name]"
echo ""
echo "See README.md for detailed usage instructions." 