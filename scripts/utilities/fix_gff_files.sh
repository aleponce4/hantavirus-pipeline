#!/bin/bash

# Script to fix GFF files for compatibility with iVar

# Process all GFF files in reference directories
for segment in L_segment M_segment S_segment; do
    segment_dir="data/references/$segment"
    
    if [ ! -d "$segment_dir" ]; then
        echo "Directory $segment_dir does not exist, skipping..."
        continue
    fi
    
    echo "Processing GFF files in $segment_dir..."
    
    # Find all GFF files
    for gff_file in "$segment_dir"/*.gff3; do
        if [ ! -f "$gff_file" ]; then
            echo "  No GFF files found in $segment_dir"
            continue
        fi
        
        echo "  Fixing $gff_file..."
        
        # Make a backup of the original file
        backup_file="${gff_file}.backup"
        cp "$gff_file" "$backup_file"
        
        # Fix line endings and any other issues
        # 1. Convert potential DOS/Windows line endings to Unix
        # 2. Ensure no wrapped lines (sometimes GFF files have wrapped lines which confuse parsers)
        
        # First, unwrap any wrapped lines and convert DOS to Unix endings
        tr -d '\r' < "$backup_file" > "$gff_file.tmp1"
        
        # Create a simpler version of the GFF file with proper tabs
        # This handles the most common issues with GFF files
        cat "$gff_file.tmp1" | awk '{
            # If it starts with # it is a comment line
            if ($0 ~ /^#/) {
                print $0;
            } 
            # If it is a data line
            else {
                # Replace spaces with tabs in the required fields
                gsub(" ","\t",$0);
                print $0;
            }
        }' > "$gff_file.tmp2"
        
        # Move fixed file into place
        mv "$gff_file.tmp2" "$gff_file"
        
        # Clean up
        rm -f "$gff_file.tmp1"
        
        # Check if the file was fixed properly
        if grep -q "^JN" "$gff_file" || grep -q "^OR" "$gff_file"; then
            echo "  GFF file fixed successfully."
        else
            echo "  WARNING: GFF file may still have issues."
            # Restore from backup if the fix seems to have broken the file
            if [ ! -s "$gff_file" ]; then
                echo "  Error: Fixed file is empty. Restoring from backup..."
                cp "$backup_file" "$gff_file"
            fi
        fi
    done
done

echo "All GFF files processed." 