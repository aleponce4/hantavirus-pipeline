#!/bin/bash
# Simplified consensus generation script using iVar

sample=$1

# Source the configuration file
source ./config.sh

if [ "$FIRST_PASS" = "true" ]; then
    echo "Generating consensus for sample: $sample (FIRST PASS)"
    RESULTS_DIR="results/first_pass"
else
    echo "Generating consensus for sample: $sample (SECOND PASS)"
    RESULTS_DIR="results/second_pass"
fi

echo "Using minimum coverage threshold: $MIN_CONSENSUS_COVERAGE"

# Process each segment
for segment in M_segment S_segment; do  # Skip L_segment - no data
    segment_dir="data/references/$segment"
    
    if [ -d "$segment_dir" ] && [ "$(ls -A $segment_dir)" ]; then
        # Use original reference for both passes (simplified)
        reference=$(ls "$segment_dir"/Reference_*.fasta | head -n 1)
        echo "Using reference for $segment: $reference"
        
        # Set up paths
        bam_file="$RESULTS_DIR/$sample/alignment/$segment/${sample}.bam"
        out_dir="$RESULTS_DIR/$sample/consensus/$segment"
        mkdir -p "$out_dir"
        
        echo "Generating consensus for $sample - $segment"
        
        # Check if we have a cleaned BAM from variant calling (after primer mismatch removal)
        cleaned_bam="$RESULTS_DIR/$sample/variants/$segment/cleaned.sorted.bam"
        if [ -f "$cleaned_bam" ]; then
            echo "  Using cleaned BAM (post primer-mismatch removal)"
            input_bam="$cleaned_bam"
        else
            # Check for primer-trimmed BAM
            trimmed_bam="$RESULTS_DIR/$sample/variants/$segment/trimmed.sorted.bam"
            if [ -f "$trimmed_bam" ]; then
                echo "  Using primer-trimmed BAM"
                input_bam="$trimmed_bam"
            else
                echo "  Using original alignment BAM"
                input_bam="$bam_file"
            fi
        fi
        
        # Generate majority consensus with iVar
        echo "  Creating majority consensus with iVar..."
        samtools mpileup -aa -A -d 0 -Q 0 --reference "$reference" "$input_bam" | \
            ivar consensus -p "$out_dir/${sample}_consensus" \
            -m "$MIN_CONSENSUS_COVERAGE" -t 0.5 -q "$MIN_VARIANT_QUALITY" -n N \
            -i "${sample}_${segment}_consensus"
        
        # Check if consensus was generated
        if [ -f "$out_dir/${sample}_consensus.fa" ]; then
            # Count coverage statistics
            total_positions=$(samtools view -H "$reference" | grep "^@SQ" | awk '{print $3}' | sed 's/LN://')
            n_count=$(grep -o "N" "$out_dir/${sample}_consensus.fa" | wc -l)
            coverage_positions=$((total_positions - n_count))
            coverage_pct=$(awk -v cov="$coverage_positions" -v total="$total_positions" 'BEGIN {printf "%.1f", cov*100/total}')
            
            echo "  Consensus generated successfully"
            echo "    Total positions: $total_positions"
            echo "    Covered positions: $coverage_positions (${coverage_pct}%)"
            echo "    Masked positions (Ns): $n_count"
            
            # Rename to expected filename for compatibility
            mv "$out_dir/${sample}_consensus.fa" "$out_dir/${sample}_consensus.fasta"
            
            echo "  Consensus saved to: $out_dir/${sample}_consensus.fasta"
        else
            echo "  ERROR: Failed to generate consensus for $sample - $segment"
            exit 1
        fi
    fi
done

# Create placeholder for L_segment
L_dir="$RESULTS_DIR/$sample/consensus/L_segment"
mkdir -p "$L_dir"
echo ">$sample L_segment consensus (NO DATA AVAILABLE)" > "$L_dir/${sample}_consensus.fasta"
echo "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" >> "$L_dir/${sample}_consensus.fasta"
echo "Created placeholder consensus for L_segment"

echo "Consensus generation completed for sample: $sample" 