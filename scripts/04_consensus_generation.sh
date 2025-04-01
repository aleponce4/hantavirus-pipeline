#!/bin/bash

# Consensus generation script

sample=$1

# Check if this is the first pass
if [ "$FIRST_PASS" = "true" ]; then
    echo "Generating consensus for sample: $sample (FIRST PASS)"
    RESULTS_DIR="results/first_pass"
    # For first pass, we don't need to compare with previous pass
    COMPARE_WITH_PREVIOUS=false
else
    echo "Generating consensus for sample: $sample (SECOND PASS)"
    RESULTS_DIR="results/second_pass"
    # For second pass, track the first pass consensus for comparison
    FIRST_PASS_DIR="results/first_pass"
    COMPARE_WITH_PREVIOUS=true
fi

# Set minimum coverage threshold - positions with coverage below this will be masked with N
MIN_COVERAGE=10
echo "Using minimum coverage threshold: $MIN_COVERAGE (positions below this will be masked with N)"

# Process each segment
for segment in L_segment M_segment S_segment; do
    segment_dir="data/references/$segment"
    
    # Check if there are any files in the segment directory
    if [ -d "$segment_dir" ] && [ "$(ls -A $segment_dir)" ]; then
        # Reference selection strategy
        if [ "$FIRST_PASS" = "true" ]; then
            # For first pass, always use the original reference
            reference=$(ls $segment_dir/*.fasta | grep -v "metaconsensus.fasta" | head -n 1)
            echo "FIRST PASS: Using original reference for $segment: $reference"
        else
            # For second pass, always use the metaconsensus reference
            metaconsensus="data/references/$segment/metaconsensus.fasta"
            if [ -f "$metaconsensus" ]; then
                reference=$metaconsensus
                echo "SECOND PASS: Using metaconsensus reference for $segment: $reference"
            else
                echo "ERROR: Metaconsensus reference not found for $segment. Cannot proceed with second pass."
                exit 1
            fi
        fi
        
        if [ -n "$reference" ]; then
            # Input files
            bam_file="$RESULTS_DIR/$sample/alignment/$segment/${sample}.bam"
            vcf_file="$RESULTS_DIR/$sample/variants/$segment/${sample}.vcf.gz"
            
            # Output directory
            out_dir="$RESULTS_DIR/$sample/consensus/$segment"
            mkdir -p $out_dir
            
            # Generate consensus FASTA
            consensus_file="$out_dir/${sample}_consensus.fasta"
            
            # Debug: Check if VCF has variants
            echo "Checking VCF file for variants..."
            variant_count=$(zcat $vcf_file | grep -v "^#" | wc -l)
            echo "Found $variant_count variants in VCF file"
            
            # Generate mask for low coverage regions
            echo "Generating coverage mask for regions below ${MIN_COVERAGE}x coverage..."
            mask_bed="$out_dir/${sample}_mask.bed"
            samtools depth -a "$bam_file" | awk -v min="$MIN_COVERAGE" '$3 < min {print $1"\t"$2-1"\t"$2}' > "$mask_bed"
            
            # Count masked positions
            mask_count=$(wc -l < "$mask_bed")
            echo "Number of positions to be masked due to low coverage: $mask_count"
            
            # Make sure reference is indexed
            if [ ! -f "${reference}.fai" ]; then
                samtools faidx "$reference"
            fi
            
            # Generate consensus with bcftools
            echo "Generating consensus using bcftools consensus..."
            temp_consensus="$out_dir/temp_consensus.fa"
            
            # Use bcftools with -I flag to handle empty IDs in LoFreq's VCF
            bcftools consensus -I -f "$reference" "$vcf_file" -o "$temp_consensus" 2> "$out_dir/bcftools.log"
            
            # Check if consensus generation succeeded
            if [ $? -ne 0 ] || [ ! -s "$temp_consensus" ]; then
                echo "ERROR: bcftools consensus failed. Check the log file for details:"
                cat "$out_dir/bcftools.log"
                echo "Failed to generate consensus for $sample - $segment"
                exit 1
            fi
            
            # Apply masking for low coverage regions
            if [ "$mask_count" -gt 0 ]; then
                echo "Applying coverage mask to consensus..."
                bedtools maskfasta -fi "$temp_consensus" -bed "$mask_bed" -fo "$consensus_file" -mc N
                
                # Check if masking succeeded
                if [ $? -ne 0 ] || [ ! -s "$consensus_file" ]; then
                    echo "WARNING: bedtools maskfasta failed. Using unmasked consensus."
                    cp "$temp_consensus" "$consensus_file"
                fi
            else
                # No masking needed
                cp "$temp_consensus" "$consensus_file"
            fi
            
            # Fix the header
            sed -i "1s/.*/>$sample $segment consensus/" "$consensus_file"
            
            # Count mismatches from the log
            mismatch_count=$(grep -c "Reference FASTA/BCF inconsistency" "$out_dir/bcftools.log")
            if [ "$mismatch_count" -gt 0 ]; then
                echo "WARNING: $mismatch_count reference mismatches detected. See $out_dir/bcftools.log for details."
                if [ "$mismatch_count" -gt 20 ]; then
                    echo "Consider checking that the reference FASTA exactly matches what was used during alignment."
                fi
            fi
            
            # Clean up temporary files
            rm -f "$temp_consensus"
            
            # For second pass, compare with first pass consensus to evaluate improvement
            if [ "$COMPARE_WITH_PREVIOUS" == "true" ]; then
                first_pass_consensus="$FIRST_PASS_DIR/$sample/consensus/$segment/${sample}_consensus.fasta"
                if [ -f "$first_pass_consensus" ]; then
                    echo "Comparing first pass and second pass consensus sequences..."
                    
                    # Create a comparison directory if it doesn't exist
                    comparison_dir="$RESULTS_DIR/$sample/comparison/$segment"
                    mkdir -p "$comparison_dir"
                    
                    # Get sequences without headers for comparison
                    grep -v "^>" "$first_pass_consensus" > "$comparison_dir/first_pass_seq.txt"
                    grep -v "^>" "$consensus_file" > "$comparison_dir/second_pass_seq.txt"
                    
                    # Basic diff to count differences
                    diff_count=$(diff -y --suppress-common-lines "$comparison_dir/first_pass_seq.txt" "$comparison_dir/second_pass_seq.txt" | wc -l)
                    echo "CONSENSUS COMPARISON: Found $diff_count differences between first and second pass"
                    
                    # Count Ns in each consensus
                    first_pass_Ns=$(grep -o "N" "$comparison_dir/first_pass_seq.txt" | wc -l)
                    second_pass_Ns=$(grep -o "N" "$comparison_dir/second_pass_seq.txt" | wc -l)
                    echo "COVERAGE COMPARISON: First pass had $first_pass_Ns Ns, Second pass has $second_pass_Ns Ns"
                    
                    # Calculate N reduction percentage
                    if [ "$first_pass_Ns" -gt 0 ]; then
                        n_reduction=$(( (first_pass_Ns - second_pass_Ns) * 100 / first_pass_Ns ))
                        echo "COVERAGE IMPROVEMENT: Second pass reduced Ns by $n_reduction% compared to first pass"
                    elif [ "$first_pass_Ns" -eq 0 ] && [ "$second_pass_Ns" -eq 0 ]; then
                        echo "COVERAGE COMPARISON: Both passes have complete coverage (no Ns)"
                    fi
                    
                    # Calculate consensus similarity
                    first_pass_length=$(wc -c < "$comparison_dir/first_pass_seq.txt")
                    similarity_pct=$(( (first_pass_length - diff_count) * 100 / first_pass_length ))
                    echo "CONSENSUS SIMILARITY: First and second pass are $similarity_pct% identical"
                    
                    # Save comparison summary
                    {
                        echo "CONSENSUS COMPARISON SUMMARY FOR $sample - $segment"
                        echo "----------------------------------------------"
                        echo "Differences between passes: $diff_count positions"
                        echo "First pass Ns: $first_pass_Ns"
                        echo "Second pass Ns: $second_pass_Ns"
                        if [ "$first_pass_Ns" -gt 0 ]; then
                            echo "N reduction: $n_reduction%"
                        fi
                        echo "Sequence similarity: $similarity_pct%"
                        echo "----------------------------------------------"
                        echo "First pass variants: $(zcat "$FIRST_PASS_DIR/$sample/variants/$segment/${sample}.vcf.gz" | grep -v "^#" | wc -l)"
                        echo "Second pass variants: $variant_count"
                        echo "----------------------------------------------"
                    } > "$comparison_dir/comparison_summary.txt"
                    
                    # Generate a more detailed diff file for inspection
                    diff -y "$comparison_dir/first_pass_seq.txt" "$comparison_dir/second_pass_seq.txt" > "$comparison_dir/full_diff.txt"
                    
                    echo "Comparison details saved to $comparison_dir/comparison_summary.txt"
                else
                    echo "WARNING: Could not find first pass consensus at $first_pass_consensus for comparison"
                fi
            fi
            
            echo "Successfully generated consensus for $sample - $segment"
        fi
    fi
done

echo "Consensus generation completed for sample: $sample" 