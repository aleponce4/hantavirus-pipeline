#!/bin/bash

# Reference refinement script - generates metaconsensus from merged variants

threads=$1

# Check if this is the first pass
if [ "$FIRST_PASS" = "true" ]; then
    echo "Reference refinement - generating metaconsensus sequences from merged variants (from FIRST PASS results)"
    RESULTS_DIR="results/first_pass"
else
    echo "Reference refinement - generating metaconsensus sequences from merged variants (from SECOND PASS results)"
    RESULTS_DIR="results"
fi

# Process each segment
for segment in L_segment M_segment S_segment; do
    segment_dir="data/references/$segment"
    
    # Check if there are any files in the segment directory
    if [ -d "$segment_dir" ] && [ "$(ls -A $segment_dir)" ]; then
        # Get the original reference
        original_reference=$(ls $segment_dir/*.fasta | grep -v "metaconsensus.fasta" | head -n 1)
        echo "Original reference for $segment: $original_reference"
        
        # Temporary directory for processing
        tmp_dir=$(mktemp -d)
        
        # Find all VCF files for this segment
        echo "Collecting variants from all samples for $segment"
        vcf_list=""
        sample_count=0
        for sample_dir in $RESULTS_DIR/*/variants/$segment; do
            if [ -d "$sample_dir" ]; then
                for vcf in $sample_dir/*.vcf.gz; do
                    if [ -f "$vcf" ]; then
                        # Make sure VCF is indexed
                        if [ ! -f "${vcf}.tbi" ]; then
                            tabix -p vcf "$vcf"
                        fi
                        vcf_list="$vcf_list $vcf"
                        ((sample_count++))
                    fi
                done
            fi
        done
        
        # If we have no VCF files, use original reference
        if [ -z "$vcf_list" ] || [ "$sample_count" -eq 0 ]; then
            echo "No VCF files found for $segment, using original reference as metaconsensus"
            cat "$original_reference" | sed "1s/>.*/>metaconsensus $segment/" > "data/references/$segment/metaconsensus.fasta"
            rm -rf "$tmp_dir"
            continue
        fi
        
        echo "Found $sample_count VCF files for $segment"
        
        # Make sure reference is indexed
        if [ ! -f "${original_reference}.fai" ]; then
            samtools faidx "$original_reference"
        fi
        
        # Direct approach to process VCFs
        if [ "$sample_count" -eq 1 ]; then
            # If only one sample, use its VCF directly
            echo "Only one sample VCF found, using it directly"
            single_vcf=$(echo $vcf_list | tr -d ' ')
            
            # Check for variants in this VCF
            variant_count=$(bcftools view -H "$single_vcf" | wc -l)
            echo "VCF contains $variant_count variants"
            
            if [ "$variant_count" -eq 0 ]; then
                echo "WARNING: No variants found in VCF. Using original reference as metaconsensus."
                cat "$original_reference" | sed "1s/>.*/>metaconsensus $segment/" > "data/references/$segment/metaconsensus.fasta"
                rm -rf "$tmp_dir"
                continue
            fi
            
            # Show sample variants
            echo "Sample variants (first 5):"
            bcftools view -H "$single_vcf" | head -5
            
            # Apply variants to reference
            echo "Applying variants to reference"
            # Use the -s '-' flag to ignore sample information and just apply all ALT alleles
            bcftools consensus -s '-' -f "$original_reference" -o "$tmp_dir/metaconsensus.fasta" "$single_vcf"
            
        else
            # For multiple samples, normalize and merge VCFs
            echo "Normalizing and merging multiple VCFs"
            
            # Normalize VCFs before merging
            echo "Step 1: Normalizing individual VCFs"
            norm_vcfs=""
            for vcf in $vcf_list; do
                base_name=$(basename "$vcf" .vcf.gz)
                norm_vcf="$tmp_dir/${base_name}.norm.vcf.gz"
                echo "  Normalizing: $vcf"
                bcftools norm -f "$original_reference" -m +any "$vcf" -Oz -o "$norm_vcf"
                tabix -p vcf "$norm_vcf"
                norm_vcfs="$norm_vcfs $norm_vcf"
            done
            
            # Merge normalized VCFs
            echo "Step 2: Merging normalized VCFs"
            merged_vcf="$tmp_dir/merged.vcf.gz"
            bcftools merge --force-samples -m none -Oz -o "$merged_vcf" $norm_vcfs
            tabix -p vcf "$merged_vcf"
            
            # Check variant count in merged VCF
            variant_count=$(bcftools view -H "$merged_vcf" | wc -l)
            echo "Merged VCF contains $variant_count variants"
            
            if [ "$variant_count" -eq 0 ]; then
                echo "WARNING: No variants found in merged VCF. Using original reference as metaconsensus."
                cat "$original_reference" | sed "1s/>.*/>metaconsensus $segment/" > "data/references/$segment/metaconsensus.fasta"
                rm -rf "$tmp_dir"
                continue
            fi
            
            # Show merged variants
            echo "Merged variants (first 5):"
            bcftools view -H "$merged_vcf" | head -5
            
            # Apply variants to reference
            echo "Step 3: Applying merged variants to reference"
            # Use the -s '-' flag to ignore sample information and just apply all ALT alleles
            bcftools consensus -s '-' -f "$original_reference" -o "$tmp_dir/metaconsensus.fasta" "$merged_vcf"
        fi
        
        # Verify consensus was created
        if [ ! -s "$tmp_dir/metaconsensus.fasta" ]; then
            echo "ERROR: Failed to generate metaconsensus. Using original reference."
            cat "$original_reference" | sed "1s/>.*/>metaconsensus $segment/" > "data/references/$segment/metaconsensus.fasta"
            rm -rf "$tmp_dir"
            continue
        fi
        
        # Fix the header
        sed -i "1s/>.*/>metaconsensus $segment/" "$tmp_dir/metaconsensus.fasta"
        
        # Verify differences from reference
        echo "Confirming changes by comparing to original reference"
        ref_seq=$(grep -v "^>" "$original_reference")
        meta_seq=$(grep -v "^>" "$tmp_dir/metaconsensus.fasta")
        
        if [ "$ref_seq" = "$meta_seq" ]; then
            echo "WARNING: Metaconsensus is identical to original reference."
            echo "No variants were successfully applied."
        else
            diff_count=$(diff <(echo "$ref_seq") <(echo "$meta_seq") | grep -c "^<\|^>")
            echo "Metaconsensus differs from original reference at $diff_count positions"
        fi
        
        # Update the reference with the new metaconsensus
        echo "Saving metaconsensus reference"
        cp "$tmp_dir/metaconsensus.fasta" "data/references/$segment/metaconsensus.fasta"
        
        echo "Metaconsensus generation complete for $segment"
        
        # Clean up temporary files
        rm -rf "$tmp_dir"
    fi
done

echo "Reference refinement completed" 