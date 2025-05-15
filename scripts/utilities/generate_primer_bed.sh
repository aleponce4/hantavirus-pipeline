#!/bin/bash

# Script to generate BED files for primer masking

# Output directory for BED files
mkdir -p data/primers/bed_files

# Process each segment
for segment in "S_segment" "M_segment" "L_segment"; do
    echo "Generating BED file for $segment..."
    
    # Define output BED file
    BED_FILE="data/primers/bed_files/${segment}_primers.bed"
    
    # Get the reference sequence for this segment
    segment_dir="data/references/$segment"
    
    # Check if the segment directory exists
    if [ ! -d "$segment_dir" ]; then
        echo "Warning: Reference directory $segment_dir does not exist. Skipping."
        continue
    fi
    
    # Find reference file
    reference_files=($(ls "$segment_dir"/*.fasta 2>/dev/null | grep -v "metaconsensus.fasta"))
    
    if [ ${#reference_files[@]} -eq 0 ]; then
        echo "Warning: No reference FASTA file found in $segment_dir. Skipping."
        continue
    fi
    
    reference="${reference_files[0]}"
    echo "Using reference: $reference"
    
    # Extract reference sequence name
    ref_name=$(grep ">" "$reference" | head -1 | sed 's/>//' | cut -d' ' -f1)
    
    if [ -z "$ref_name" ]; then
        echo "Warning: Could not extract reference name from $reference. Using segment name instead."
        ref_name="${segment}"
    fi
    
    echo "Reference sequence name: $ref_name"
    
    # Clear previous BED file if it exists
    > "$BED_FILE"
    
    # Parse the primers file
    while IFS=, read -r fname fregion fseq flen fgc ftm region rname rregion rseq rlen rgc rtm rest; do
        # Skip header line
        if [[ $fname == "F Name" ]]; then
            continue
        fi
        
        # Process only primers for the current segment
        if [[ ($fname == S* && $segment == "S_segment") || 
              ($fname == M* && $segment == "M_segment") || 
              ($fname == L* && $segment == "L_segment") ]]; then
            
            # Extract positions from region names
            # Format is typically like "1-32F" for forward and "452-480R" for reverse
            
            # Forward primer info
            f_start=$(echo "$fregion" | cut -d'-' -f1)
            f_end=$(echo "$fregion" | cut -d'-' -f2 | sed 's/[^0-9]//g')
            
            # Adjust to 0-based BED format (BED is 0-based, primers are 1-based)
            f_start=$((f_start - 1))
            
            # Reverse primer info
            r_start=$(echo "$rregion" | cut -d'-' -f1)
            r_end=$(echo "$rregion" | cut -d'-' -f2 | sed 's/[^0-9]//g')
            
            # Adjust to 0-based BED format
            r_start=$((r_start - 1))
            
            # Write primers to BED file
            echo -e "$ref_name\t$f_start\t$f_end\t${fname}_primer" >> "$BED_FILE"
            echo -e "$ref_name\t$r_start\t$r_end\t${rname}_primer" >> "$BED_FILE"
            
            echo "  Added primer pair: $fname / $rname"
        fi
    done < "data/primers/Primers.csv"
    
    # Count entries in BED file
    primer_count=$(wc -l < "$BED_FILE")
    echo "Created BED file for $segment with $primer_count primer regions"
done

echo "BED files generated in data/primers/bed_files/" 