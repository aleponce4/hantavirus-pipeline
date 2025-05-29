#!/bin/bash
# Consensus generation script - uses consensus from full iVar pipeline

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
        # Use original reference
        reference=$(ls "$segment_dir"/Reference_*.fasta | head -n 1)
        echo "Using reference for $segment: $reference"
        
        # Set up paths
        variants_dir="$RESULTS_DIR/$sample/variants/$segment"
        out_dir="$RESULTS_DIR/$sample/consensus/$segment"
        mkdir -p "$out_dir"
        
        echo "Generating consensus for $sample - $segment"
        
        # Check if consensus was already generated during variant calling
        if [ -f "$variants_dir/consensus.fa" ]; then
            echo "  Using consensus generated during variant calling pipeline"
            cp "$variants_dir/consensus.fa" "$out_dir/${sample}_consensus.fasta"
            
            # Fix header to match expected format
            sed -i "1s/.*/>$sample $segment consensus/" "$out_dir/${sample}_consensus.fasta"
            
        else
            echo "  Generating new consensus from final alignment..."
            
            # Find the best available BAM file (cleaned > final > trimmed > original)
            if [ -f "$variants_dir/final_alignment.bam" ]; then
                input_bam="$variants_dir/final_alignment.bam"
                echo "  Using final cleaned alignment"
            elif [ -f "$variants_dir/cleaned.sorted.bam" ]; then
                input_bam="$variants_dir/cleaned.sorted.bam"
                echo "  Using cleaned BAM"
            elif [ -f "$variants_dir/trimmed.sorted.bam" ]; then
                input_bam="$variants_dir/trimmed.sorted.bam"
                echo "  Using primer-trimmed BAM"
            else
                bam_file="$RESULTS_DIR/$sample/alignment/$segment/${sample}.bam"
                input_bam="$bam_file"
                echo "  Using original alignment BAM"
            fi
            
            # Generate consensus with iVar
            samtools mpileup -aa -A -d 0 -Q 0 --reference "$reference" "$input_bam" | \
                ivar consensus -p "$out_dir/${sample}_consensus" \
                -m "$MIN_CONSENSUS_COVERAGE" -t 0.5 -q "$MIN_VARIANT_QUALITY" -n N \
                -i "${sample}_${segment}_consensus"
            
            # Rename to expected filename
            if [ -f "$out_dir/${sample}_consensus.fa" ]; then
                mv "$out_dir/${sample}_consensus.fa" "$out_dir/${sample}_consensus.fasta"
            fi
        fi
        
        # Check if consensus was generated successfully
        if [ -f "$out_dir/${sample}_consensus.fasta" ]; then
            # Count coverage statistics
            total_positions=$(samtools view -H "$reference" | grep "^@SQ" | awk '{print $3}' | sed 's/LN://')
            n_count=$(grep -o "N" "$out_dir/${sample}_consensus.fasta" | wc -l)
            coverage_positions=$((total_positions - n_count))
            
            # Fix division by zero
            if [ "$total_positions" -gt 0 ]; then
                coverage_pct=$(awk -v cov="$coverage_positions" -v total="$total_positions" 'BEGIN {printf "%.1f", cov*100/total}')
            else
                coverage_pct="0.0"
            fi
            
            echo "  Consensus generated successfully"
            echo "    Total positions: $total_positions"
            echo "    Covered positions: $coverage_positions (${coverage_pct}%)"
            echo "    Masked positions (Ns): $n_count"
            
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