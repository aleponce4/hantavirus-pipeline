#!/bin/bash

# Primer evaluation script for Hantavirus pipeline
# Evaluates the quality of primers in Primers.csv against various criteria

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

# Generate improvement suggestions for the S segment primers
echo -e "\nGenerating improvement suggestions for S segment primers..." | tee -a $PRIMER_EVAL_LOG
python scripts/primer_evaluation/improve_primers.py \
    --eval_report "$RESULTS_DIR/S_segment_primer_evaluation.csv" \
    --reference "$CONSENSUS" \
    --output "$RESULTS_DIR/S_segment_primer_improvements.csv" \
    --blast_db "$BLAST_DB" \
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

# Generate improvement suggestions for the M segment primers
echo -e "\nGenerating improvement suggestions for M segment primers..." | tee -a $PRIMER_EVAL_LOG
python scripts/primer_evaluation/improve_primers.py \
    --eval_report "$RESULTS_DIR/M_segment_primer_evaluation.csv" \
    --reference "$CONSENSUS" \
    --output "$RESULTS_DIR/M_segment_primer_improvements.csv" \
    --blast_db "$BLAST_DB" \
    2>&1 | tee -a $PRIMER_EVAL_LOG

# Add code for L segment when available

echo -e "\nPrimer evaluation completed at $(date)" | tee -a $PRIMER_EVAL_LOG
echo "Results available in $RESULTS_DIR"

# Summarize improvement suggestions
echo -e "\n=== PRIMER IMPROVEMENT SUMMARY ===" | tee -a $PRIMER_EVAL_LOG

# S segment summary
if [ -f "$RESULTS_DIR/S_segment_primer_improvements.csv" ]; then
    IMPROVED_COUNT=$(wc -l < "$RESULTS_DIR/S_segment_primer_improvements.csv")
    # Subtract 1 for the header line
    if [ $IMPROVED_COUNT -gt 1 ]; then
        IMPROVED_COUNT=$((IMPROVED_COUNT - 1))
        echo "S segment: Generated improvements for $IMPROVED_COUNT primers." | tee -a $PRIMER_EVAL_LOG
        echo "Full details in: $RESULTS_DIR/S_segment_primer_improvements.csv" | tee -a $PRIMER_EVAL_LOG
        
        # Generate a quick summary table of changed primers
        if command -v column > /dev/null; then
            echo -e "\nSummary of changed S segment primers:" | tee -a $PRIMER_EVAL_LOG
            (echo "Primer,Type,Changes,Length Δ"; tail -n +2 "$RESULTS_DIR/S_segment_primer_improvements.csv" | cut -d, -f1,2,15,16 | sort -t, -k3,3nr) | column -t -s, | tee -a $PRIMER_EVAL_LOG
        fi
    else
        echo "S segment: No primer improvements needed - all primers pass quality checks." | tee -a $PRIMER_EVAL_LOG
    fi
fi

# M segment summary
if [ -f "$RESULTS_DIR/M_segment_primer_improvements.csv" ]; then
    IMPROVED_COUNT=$(wc -l < "$RESULTS_DIR/M_segment_primer_improvements.csv")
    # Subtract 1 for the header line
    if [ $IMPROVED_COUNT -gt 1 ]; then
        IMPROVED_COUNT=$((IMPROVED_COUNT - 1))
        echo -e "\nM segment: Generated improvements for $IMPROVED_COUNT primers." | tee -a $PRIMER_EVAL_LOG
        echo "Full details in: $RESULTS_DIR/M_segment_primer_improvements.csv" | tee -a $PRIMER_EVAL_LOG
        
        # Generate a quick summary table of changed primers
        if command -v column > /dev/null; then
            echo -e "\nSummary of changed M segment primers:" | tee -a $PRIMER_EVAL_LOG
            (echo "Primer,Type,Changes,Length Δ"; tail -n +2 "$RESULTS_DIR/M_segment_primer_improvements.csv" | cut -d, -f1,2,15,16 | sort -t, -k3,3nr) | column -t -s, | tee -a $PRIMER_EVAL_LOG
        fi
    else
        echo "M segment: No primer improvements needed - all primers pass quality checks." | tee -a $PRIMER_EVAL_LOG
    fi
fi 