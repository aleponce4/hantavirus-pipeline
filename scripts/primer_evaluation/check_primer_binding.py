#!/usr/bin/env python3
"""
Script to check if primers bind to their intended sites on a reference sequence.
This standalone script can be used to quickly check primer binding without running
the full evaluation pipeline.
"""

import argparse
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import pandas as pd

def verify_primer_binding_site(primer_seq, reference_seq, position=None, is_reverse=False, flanking_region=40):
    """
    Verify if a primer matches its intended binding region in the reference sequence.
    
    Args:
        primer_seq (str): Primer sequence
        reference_seq (str): Full reference sequence
        position (tuple): Tuple containing (start, end) positions of primer, or None to search whole sequence
        is_reverse (bool): Whether this is a reverse primer
        flanking_region (int): Number of bases to check on either side of the expected position
    
    Returns:
        dict: Dictionary with match percentage, actual position, and match status
    """
    # For reverse primers, we need to use the reverse complement for matching
    primer_to_match = primer_seq
    if is_reverse:
        primer_to_match = str(Seq(primer_seq).reverse_complement())
    
    # If position is provided, check only that region (with flanking)
    if position is not None:
        start, end = position
        
        # Adjust for 0-based indexing
        start = max(0, start)
        end = min(len(reference_seq) - 1, end)
        
        # Expand region to check for slight positional variations
        check_start = max(0, start - flanking_region)
        check_end = min(len(reference_seq), end + flanking_region + len(primer_seq))
        
        # Extract the region to check
        region_to_check = reference_seq[check_start:check_end]
    else:
        # Check the entire sequence
        region_to_check = reference_seq
        check_start = 0
    
    # Find best match within the region
    best_match = 0
    best_pos = None
    
    # Slide the primer across the region
    for i in range(len(region_to_check) - len(primer_to_match) + 1):
        match_count = sum(a == b for a, b in zip(primer_to_match, region_to_check[i:i+len(primer_to_match)]))
        match_percent = match_count / len(primer_to_match) * 100
        
        if match_percent > best_match:
            best_match = match_percent
            best_pos = check_start + i
    
    # Determine if it's a good match (over 80% is generally acceptable for PCR)
    match_status = "Good" if best_match >= 80 else "Poor"
    
    result = {
        "match_percent": best_match,
        "actual_position": (best_pos, best_pos + len(primer_to_match) - 1) if best_pos is not None else None,
        "status": match_status
    }
    
    # Add position offset information if we had an expected position
    if position is not None and best_pos is not None:
        result["distance_from_expected"] = best_pos - start
    
    return result

def run_blast_and_get_evalue(primer_seq, blast_db):
    """Run BLAST on a primer sequence and get the best E-value"""
    # Create temporary file
    temp_file = "temp_primer.fasta"
    with open(temp_file, 'w') as f:
        f.write(f">primer\n{primer_seq}\n")
    
    # Run BLAST
    cmd = [
        "blastn", 
        "-query", temp_file, 
        "-db", blast_db, 
        "-outfmt", "6", 
        "-evalue", "10", 
        "-word_size", "7", 
        "-task", "blastn-short", 
        "-dust", "no"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Parse output
        if result.stdout.strip():
            # Get the best hit (first line)
            best_hit = result.stdout.strip().split('\n')[0]
            fields = best_hit.split('\t')
            e_value = float(fields[10])
            return e_value
        else:
            return float('inf')  # No hits
    except Exception as e:
        print(f"Error running BLAST: {e}")
        return float('inf')
    finally:
        # Clean up
        if os.path.exists(temp_file):
            os.remove(temp_file)

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description='Check primer binding to reference sequence')
    parser.add_argument('--primer', required=True, help='Primer sequence')
    parser.add_argument('--reference', required=True, help='Path to reference FASTA file')
    parser.add_argument('--position', help='Expected position (format: start-end)')
    parser.add_argument('--reverse', action='store_true', help='Primer is a reverse primer')
    parser.add_argument('--blast_db', help='Path to BLAST database for E-value calculation')
    
    args = parser.parse_args()
    
    # Load reference sequence
    try:
        reference_seq = str(SeqIO.read(args.reference, "fasta").seq)
    except Exception as e:
        print(f"Error reading reference file: {e}")
        sys.exit(1)
    
    # Parse position if provided
    position = None
    if args.position:
        try:
            parts = args.position.split('-')
            start = int(parts[0]) - 1  # Convert to 0-based
            end = int(''.join(c for c in parts[1] if c.isdigit())) - 1
            position = (start, end)
        except Exception as e:
            print(f"Error parsing position: {e}")
            print("Using position=None to search the entire sequence")
    
    # Check binding
    result = verify_primer_binding_site(
        args.primer,
        reference_seq,
        position,
        is_reverse=args.reverse,
        flanking_region=5
    )
    
    # Print results
    print("\nPrimer Binding Analysis")
    print("======================")
    print(f"Primer: {args.primer}")
    print(f"Length: {len(args.primer)} bp")
    print(f"Type: {'Reverse' if args.reverse else 'Forward'}")
    
    if position:
        print(f"Expected position: {position[0]+1}-{position[1]+1}")
    
    print(f"\nMatch percentage: {result['match_percent']:.1f}%")
    
    if result['actual_position']:
        actual_start, actual_end = result['actual_position']
        print(f"Actual position: {actual_start+1}-{actual_end+1}")
        
        if 'distance_from_expected' in result:
            print(f"Offset from expected: {result['distance_from_expected']} bp")
    
    print(f"Binding status: {result['status']}")
    
    # Run BLAST if database is provided
    if args.blast_db:
        e_value = run_blast_and_get_evalue(args.primer, args.blast_db)
        print(f"\nBLAST E-value: {e_value}")
        print(f"BLAST quality: {'Good' if e_value <= 1e-4 else 'Poor'}")
    
    # If we have a good match, show the alignment
    if result['actual_position']:
        actual_start, actual_end = result['actual_position']
        if result['match_percent'] > 50:  # Only show if decent match
            print("\nAlignment:")
            
            # Get the matching region from the reference
            target_seq = reference_seq[actual_start:actual_end+1]
            
            # For reverse primers, we compare against the reverse complement
            primer_for_comparison = args.primer
            if args.reverse:
                primer_for_comparison = str(Seq(args.primer).reverse_complement())
                
            # Build the alignment visual
            match_str = ''.join('|' if a == b else ' ' for a, b in zip(primer_for_comparison, target_seq))
            
            print(f"Primer: {primer_for_comparison}")
            print(f"Match:  {match_str}")
            print(f"Target: {target_seq}")

if __name__ == "__main__":
    main() 