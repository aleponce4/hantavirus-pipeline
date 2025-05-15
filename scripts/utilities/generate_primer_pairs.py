#!/usr/bin/env python3
"""
Generate primer pair information files required for iVar getmasked and removereads functions.
These files map the left (forward) and right (reverse) primers for each amplicon.
"""

import os
import csv
import re
import sys

def extract_primer_number(primer_name):
    """Extract the numeric part from primer names like SF123 or MR456."""
    match = re.search(r'[A-Z][FR](\d+)', primer_name)
    if match:
        return int(match.group(1))
    return None

def generate_primer_pairs(primers_csv, output_dir):
    """Generate primer pair information files for each segment."""
    # Dictionary to store primers by segment
    primers_by_segment = {'S': [], 'M': [], 'L': []}
    
    # Read primer information from CSV
    with open(primers_csv, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)  # Skip header
        
        for row in reader:
            if len(row) < 3:  # Skip empty rows
                continue
                
            f_name = row[0].strip()
            r_name = row[1].strip() if len(row) > 1 else ""
            
            # Determine segment based on primer name prefix
            segment = None
            if f_name.startswith('S'):
                segment = 'S'
            elif f_name.startswith('M'):
                segment = 'M'
            elif f_name.startswith('L'):
                segment = 'L'
                
            if segment and f_name and r_name:
                primers_by_segment[segment].append((f_name, r_name))
    
    # Sort primer pairs by number for each segment
    for segment, pairs in primers_by_segment.items():
        sorted_pairs = sorted(pairs, key=lambda x: extract_primer_number(x[0]))
        
        # Write to output file
        output_file = os.path.join(output_dir, f"{segment}_segment_primer_pairs.tsv")
        with open(output_file, 'w') as f:
            for forward, reverse in sorted_pairs:
                f.write(f"{forward}\t{reverse}\n")
        
        print(f"Created {segment} segment primer pair file: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <primers.csv> <output_directory>")
        sys.exit(1)
    
    primers_csv = sys.argv[1]
    output_dir = sys.argv[2]
    
    if not os.path.exists(primers_csv):
        print(f"Error: Primers CSV file '{primers_csv}' not found")
        sys.exit(1)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    generate_primer_pairs(primers_csv, output_dir)
    print("Done generating primer pair information files.") 