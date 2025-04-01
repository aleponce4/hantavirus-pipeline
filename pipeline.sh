#!/bin/bash

# Main pipeline script for viral sequencing data processing

# Create logs directory
mkdir -p logs

# Define module-specific log files
PREPROCESSING_LOG="logs/01_preprocessing.log"
ALIGNMENT_LOG="logs/02_alignment.log"
VARIANT_CALLING_LOG="logs/03_variant_calling.log"
CONSENSUS_LOG="logs/04_consensus.log"
REFINEMENT_LOG="logs/05_refinement.log"
SECOND_PASS_LOG="logs/06_second_pass.log"
PLOTS_LOG="logs/07_plots.log"
NEGATIVE_SAMPLES_LOG="logs/negative_samples.log"

# Clear any existing log files
> $PREPROCESSING_LOG
> $ALIGNMENT_LOG
> $VARIANT_CALLING_LOG
> $CONSENSUS_LOG
> $REFINEMENT_LOG
> $SECOND_PASS_LOG
> $PLOTS_LOG
> $NEGATIVE_SAMPLES_LOG

# Ensure conda environment is properly activated
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate De_Novo_pipeline

# Verify conda environment is active
if [[ "${CONDA_DEFAULT_ENV}" != "De_Novo_pipeline" ]]; then
    echo "Error: Failed to activate De_Novo_pipeline conda environment"
    exit 1
fi
echo "Successfully activated conda environment: ${CONDA_DEFAULT_ENV}"
echo "CONDA_PREFIX: $CONDA_PREFIX"

# Create results directory structure
mkdir -p results
mkdir -p results/first_pass
mkdir -p results/negative_samples
mkdir -p results/second_pass # New dedicated directory for second pass results

# Export the conda environment so it's available to subprocesses
export CONDA_PREFIX
export PATH="${CONDA_PREFIX}/bin:$PATH"

# Get unique sample names from filenames (remove _R1.fastq.gz or _R2.fastq.gz)
cd data/raw_reads
samples=$(ls *.fastq.gz | sed 's/_R[12]\.fastq\.gz$//' | sort | uniq)
cd ../../

# Calculate optimal thread counts based on available cores
TOTAL_THREADS=$(($(nproc) - 2))
echo "Total available CPU threads: ${TOTAL_THREADS}"

# Different tools have different optimal thread usage:
# - Trim Galore: max 8 cores recommended
# - BWA: Can use many cores effectively
# - Samtools: Can use many cores effectively  
# - LoFreq: Use 8 threads (hardcoded in script)
# - MAFFT: Can use many cores effectively

TRIM_THREADS=$(( TOTAL_THREADS > 8 ? 8 : TOTAL_THREADS ))
BWA_THREADS=$TOTAL_THREADS
SAMTOOLS_THREADS=$TOTAL_THREADS
LOFREQ_THREADS=4  # Fixed at 4
MAFFT_THREADS=$TOTAL_THREADS

export TRIM_THREADS BWA_THREADS SAMTOOLS_THREADS LOFREQ_THREADS MAFFT_THREADS

echo "Thread allocation:"
echo "  - Trim Galore: ${TRIM_THREADS} threads"
echo "  - BWA: ${BWA_THREADS} threads" 
echo "  - Samtools: ${SAMTOOLS_THREADS} threads"
echo "  - LoFreq: ${LOFREQ_THREADS} threads"
echo "  - MAFFT: ${MAFFT_THREADS} threads"

echo "===== FIRST PASS: Using original reference ====="

# -------------------------------------------------------
# FIRST PASS: Process all samples against original reference
# -------------------------------------------------------

# Remove any existing metaconsensus files to ensure first pass uses original references
echo "Removing any existing metaconsensus files to ensure first pass uses original references..."
for segment in L_segment M_segment S_segment; do
    metaconsensus="data/references/$segment/metaconsensus.fasta"
    if [ -f "$metaconsensus" ]; then
        echo "  Removing $metaconsensus"
        rm "$metaconsensus"
    fi
done

# Pipeline steps - FIRST PASS
for sample in $samples; do
    echo "Processing sample (first pass): $sample"
    
    # Create results directory for this sample (first pass)
    mkdir -p results/first_pass/$sample
    
    # Step 1: Preprocessing - Trim raw reads (we only do this once)
    # Store in a common location both passes can access
    mkdir -p results/trimmed/$sample
    echo "  Starting preprocessing for $sample"
    (echo "==== Processing sample: $sample ====" >> $PREPROCESSING_LOG)
    bash -c "source ${CONDA_BASE}/etc/profile.d/conda.sh && conda activate De_Novo_pipeline && OUT_DIR=results/trimmed/$sample bash scripts/01_preprocessing.sh $sample ${TRIM_THREADS}" >> $PREPROCESSING_LOG 2>&1
    echo "  Completed preprocessing for $sample"
    
    # Step 2: Reference alignment (first pass - using original reference)
    echo "  Starting reference alignment for $sample (first pass)"
    (echo "==== Processing sample: $sample (first pass) ====" >> $ALIGNMENT_LOG)
    bash -c "source ${CONDA_BASE}/etc/profile.d/conda.sh && conda activate De_Novo_pipeline && FIRST_PASS=true TRIMMED_DIR=results/trimmed/$sample bash scripts/02_reference_alignment.sh $sample ${BWA_THREADS} ${SAMTOOLS_THREADS}" >> $ALIGNMENT_LOG 2>&1
    echo "  Completed reference alignment for $sample (first pass)"
    
    # Step 3: Variant calling (first pass)
    echo "  Starting variant calling for $sample (first pass)"
    (echo "==== Processing sample: $sample (first pass) ====" >> $VARIANT_CALLING_LOG)
    bash -c "source ${CONDA_BASE}/etc/profile.d/conda.sh && conda activate De_Novo_pipeline && FIRST_PASS=true bash scripts/03_variant_calling.sh $sample" >> $VARIANT_CALLING_LOG 2>&1
    echo "  Completed variant calling for $sample (first pass)"
    
    # Step 4: Consensus generation (first pass)
    echo "  Starting consensus generation for $sample (first pass)"
    (echo "==== Processing sample: $sample (first pass) ====" >> $CONSENSUS_LOG)
    bash -c "source ${CONDA_BASE}/etc/profile.d/conda.sh && conda activate De_Novo_pipeline && FIRST_PASS=true bash scripts/04_consensus_generation.sh $sample" >> $CONSENSUS_LOG 2>&1
    echo "  Completed consensus generation for $sample (first pass)"
done

# Detect negative samples based on coverage
echo "Detecting negative samples..."
(echo "==== Detecting negative samples ====" >> $NEGATIVE_SAMPLES_LOG)
bash -c "source ${CONDA_BASE}/etc/profile.d/conda.sh && conda activate De_Novo_pipeline && bash scripts/detect_negative_samples.sh 50" >> $NEGATIVE_SAMPLES_LOG 2>&1
echo "Negative sample detection complete"

# Create a list of negative samples that should be skipped
NEGATIVE_SAMPLES_FILE="results/negative_samples/negative_samples.txt"
negative_samples=()
if [ -f "$NEGATIVE_SAMPLES_FILE" ]; then
    # First, display the content of the file for debugging
    echo "DEBUG: Content of negative samples file:"
    cat "$NEGATIVE_SAMPLES_FILE"
    
    # Populate the array with sample names from the file
    while IFS= read -r line; do
        # Only add non-empty lines
        if [ -n "$line" ]; then
            negative_samples+=("$line")
            echo "DEBUG: Added '$line' to negative_samples array"
        fi
    done < "$NEGATIVE_SAMPLES_FILE"
    
    if [ ${#negative_samples[@]} -gt 0 ]; then
        echo "Found ${#negative_samples[@]} negative samples that will skip second pass processing:"
        printf "  %s\n" "${negative_samples[@]}"
    else
        echo "WARNING: Negative samples file exists but no samples were read from it"
    fi
fi

# Step 5: Reference refinement (create metaconsensus from first pass)
echo "Starting reference refinement (creating metaconsensus)"
bash -c "source ${CONDA_BASE}/etc/profile.d/conda.sh && conda activate De_Novo_pipeline && FIRST_PASS=true bash scripts/05_reference_refinement.sh ${MAFFT_THREADS}" >> $REFINEMENT_LOG 2>&1
echo "Completed reference refinement (metaconsensus created in data/references/*/metaconsensus.fasta)"

# -------------------------------------------------------
# SECOND PASS: Reprocess all samples using the newly created metaconsensus
# -------------------------------------------------------
echo "===== SECOND PASS: Using metaconsensus reference ====="

# Pipeline steps - SECOND PASS
for sample in $samples; do
    # Skip negative samples in second pass - using improved array checking
    is_negative=false
    for neg in "${negative_samples[@]}"; do
        if [[ "$sample" == "$neg" ]]; then
            is_negative=true
            break
        fi
    done
    
    if [ "$is_negative" = true ]; then
        echo "Skipping second pass for negative sample: $sample"
        continue
    fi

    echo "Processing sample (second pass): $sample"
    
    # Create results directory for this sample (second pass)
    mkdir -p results/second_pass/$sample
    
    # Step 2: Reference alignment (second pass - using metaconsensus)
    echo "  Starting reference alignment for $sample (second pass)"
    (echo "==== Processing sample: $sample (second pass) ====" >> $SECOND_PASS_LOG)
    bash -c "source ${CONDA_BASE}/etc/profile.d/conda.sh && conda activate De_Novo_pipeline && TRIMMED_DIR=results/trimmed/$sample bash scripts/02_reference_alignment.sh $sample ${BWA_THREADS} ${SAMTOOLS_THREADS}" >> $SECOND_PASS_LOG 2>&1
    echo "  Completed reference alignment for $sample (second pass)"
    
    # Step 3: Variant calling (second pass)
    echo "  Starting variant calling for $sample (second pass)"
    (echo "==== Processing sample: $sample (second pass) ====" >> $SECOND_PASS_LOG)
    bash -c "source ${CONDA_BASE}/etc/profile.d/conda.sh && conda activate De_Novo_pipeline && bash scripts/03_variant_calling.sh $sample" >> $SECOND_PASS_LOG 2>&1
    echo "  Completed variant calling for $sample (second pass)"
    
    # Step 4: Consensus generation (second pass)
    echo "  Starting consensus generation for $sample (second pass)"
    (echo "==== Processing sample: $sample (second pass) ====" >> $SECOND_PASS_LOG)
    bash -c "source ${CONDA_BASE}/etc/profile.d/conda.sh && conda activate De_Novo_pipeline && bash scripts/04_consensus_generation.sh $sample" >> $SECOND_PASS_LOG 2>&1
    echo "  Completed consensus generation for $sample (second pass)"
done

# Step 6: Generate summary plots and analyses
echo "Generating summary plots and analyses"
(echo "==== Generating summary plots and analyses ====" >> $PLOTS_LOG)
bash -c "source ${CONDA_BASE}/etc/profile.d/conda.sh && conda activate De_Novo_pipeline && bash scripts/06_generate_plots.sh" >> $PLOTS_LOG 2>&1
echo "Summary plots and analyses complete. Results in results/plots/"

echo "Pipeline completed successfully!" 