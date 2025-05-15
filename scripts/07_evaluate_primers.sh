#!/bin/bash

# Primer evaluation script for Hantavirus pipeline
# Evaluates the quality of primers in Primers.csv against various criteria
# NOTE: Primer suggestion functionality has been removed as requested

# Create logs directory if it doesn't exist
mkdir -p logs

# Log file for this script
PRIMER_EVAL_LOG="logs/07_primer_evaluation.log"
> $PRIMER_EVAL_LOG

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

# Define file paths
PRIMERS_CSV="data/primers/Primers.csv"
RESULTS_DIR="results/primer_evaluation"
mkdir -p $RESULTS_DIR

# Create temp directory
mkdir -p $RESULTS_DIR/temp

echo "Starting primer evaluation at $(date)" | tee -a $PRIMER_EVAL_LOG

# Process S segment
echo "Evaluating S segment primers..." | tee -a $PRIMER_EVAL_LOG
BLAST_DB="data/blast_db/S_segment_blast_db"
CONSENSUS="data/references/S_segment/metaconsensus.fasta"

python scripts/primer_evaluation/evaluate_primers.py \
    --primers "$PRIMERS_CSV" \
    --segment S_segment \
    --blast_db "$BLAST_DB" \
    --consensus "$CONSENSUS" \
    --output "$RESULTS_DIR" \
    2>&1 | tee -a $PRIMER_EVAL_LOG

# Process M segment
echo -e "\nEvaluating M segment primers..." | tee -a $PRIMER_EVAL_LOG
BLAST_DB="data/blast_db/M_segment_blast_db"
CONSENSUS="data/references/M_segment/metaconsensus.fasta"

python scripts/primer_evaluation/evaluate_primers.py \
    --primers "$PRIMERS_CSV" \
    --segment M_segment \
    --blast_db "$BLAST_DB" \
    --consensus "$CONSENSUS" \
    --output "$RESULTS_DIR" \
    2>&1 | tee -a $PRIMER_EVAL_LOG

# Add code for L segment when available

echo -e "\nPrimer evaluation completed at $(date)" | tee -a $PRIMER_EVAL_LOG
echo "Results available in $RESULTS_DIR"

# Summarize evaluation results
echo -e "\n=== PRIMER EVALUATION SUMMARY ===" | tee -a $PRIMER_EVAL_LOG

# S segment summary
if [ -f "$RESULTS_DIR/S_segment_primer_evaluation.csv" ]; then
    TOTAL_COUNT=$(wc -l < "$RESULTS_DIR/S_segment_primer_evaluation.csv")
    # Subtract 1 for the header line
    TOTAL_COUNT=$((TOTAL_COUNT - 1))
    FAIL_COUNT=$(grep -c "FAIL" "$RESULTS_DIR/S_segment_primer_evaluation.csv")
    
    echo "S segment: Evaluated $TOTAL_COUNT primers, $FAIL_COUNT have issues." | tee -a $PRIMER_EVAL_LOG
    
    if [ $FAIL_COUNT -gt 0 ]; then
        echo "Problematic S segment primers:" | tee -a $PRIMER_EVAL_LOG
        grep "FAIL" "$RESULTS_DIR/S_segment_primer_evaluation.csv" | cut -d, -f1,2,3,4,5,6 | sed 's/,/ | /g' | tee -a $PRIMER_EVAL_LOG
    fi
fi

# M segment summary
if [ -f "$RESULTS_DIR/M_segment_primer_evaluation.csv" ]; then
    TOTAL_COUNT=$(wc -l < "$RESULTS_DIR/M_segment_primer_evaluation.csv")
    # Subtract 1 for the header line
    TOTAL_COUNT=$((TOTAL_COUNT - 1))
    FAIL_COUNT=$(grep -c "FAIL" "$RESULTS_DIR/M_segment_primer_evaluation.csv")
    
    echo -e "\nM segment: Evaluated $TOTAL_COUNT primers, $FAIL_COUNT have issues." | tee -a $PRIMER_EVAL_LOG
    
    if [ $FAIL_COUNT -gt 0 ]; then
        echo "Problematic M segment primers:" | tee -a $PRIMER_EVAL_LOG
        grep "FAIL" "$RESULTS_DIR/M_segment_primer_evaluation.csv" | cut -d, -f1,2,3,4,5,6 | sed 's/,/ | /g' | tee -a $PRIMER_EVAL_LOG
    fi
fi 