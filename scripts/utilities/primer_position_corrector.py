#!/usr/bin/env python3
"""
Primer Position Correction Script

This script checks the primer binding positions against the metaconsensus sequences
and creates a corrected version of the primers CSV file with updated positions.
"""

import os
import sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import re

# Define the path to the project directory
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Check and correct primer binding positions.')
    parser.add_argument('--primers', required=True, help='Path to primers CSV file')
    parser.add_argument('--output', required=True, help='Path to output corrected CSV file')
    parser.add_argument('--s_consensus', help='Path to S segment metaconsensus file', 
                      default=os.path.join(BASE_DIR, 'data/references/S_segment/metaconsensus.fasta'))
    parser.add_argument('--m_consensus', help='Path to M segment metaconsensus file',
                      default=os.path.join(BASE_DIR, 'data/references/M_segment/metaconsensus.fasta'))
    parser.add_argument('--min_match', type=float, default=0.8, help='Minimum match percentage (0.0-1.0)')
    parser.add_argument('--min_diff', type=int, default=5, help='Minimum position difference to correct')
    return parser.parse_args()

def load_primers(primers_file):
    """Load primers from CSV file."""
    print(f"Loading primers from {primers_file}")
    try:
        df = pd.read_csv(primers_file, encoding='utf-8')
        
        # Clean column names
        df.columns = df.columns.str.strip()
        clean_col_names = {}
        for col in df.columns:
            clean_name = col.replace('´', "'").replace('°C', ' C').replace('ºC', ' C')
            if clean_name != col:
                clean_col_names[col] = clean_name
        
        if clean_col_names:
            df = df.rename(columns=clean_col_names)
        
        return df
    except Exception as e:
        print(f"Error loading primers: {e}")
        sys.exit(1)

def load_consensus_sequences(s_path, m_path):
    """Load consensus sequences from fasta files."""
    try:
        s_seq = str(SeqIO.read(s_path, "fasta").seq)
        m_seq = str(SeqIO.read(m_path, "fasta").seq)
        return {'S': s_seq, 'M': m_seq}
    except Exception as e:
        print(f"Error loading consensus sequences: {e}")
        sys.exit(1)

def verify_primer_binding(primer_seq, consensus_seq, is_reverse=False, min_match=0.8):
    """Find the actual binding site of a primer in the reference sequence."""
    # Get reverse complement if it's a reverse primer
    if is_reverse:
        search_seq = str(Seq(primer_seq).reverse_complement())
    else:
        search_seq = primer_seq
    
    # Clean the sequence (remove spaces, etc.)
    search_seq = ''.join(c for c in search_seq if c.upper() in 'ACGT')
    
    # Find all close matches (allowing for some mismatches)
    best_match_pos = -1
    best_match_score = 0
    
    # Slide the primer across the entire sequence
    for i in range(len(consensus_seq) - len(search_seq) + 1):
        ref_segment = consensus_seq[i:i+len(search_seq)]
        
        # Count matches
        match_count = sum(a == b for a, b in zip(search_seq, ref_segment))
        match_score = match_count / len(search_seq)
        
        if match_score > best_match_score:
            best_match_score = match_score
            best_match_pos = i
    
    if best_match_pos >= 0 and best_match_score >= min_match:
        # Convert to 1-based positions for the CSV file
        start_pos = best_match_pos + 1
        end_pos = start_pos + len(search_seq) - 1
        
        return (start_pos, end_pos, best_match_score)
    else:
        return None

def update_region_format(region_str, start, end, is_reverse=False):
    """Update the region string with new positions."""
    if is_reverse:
        return f"{start}-{end}R"
    else:
        return f"{start}-{end}F"

def parse_region(region_str):
    """Parse a region string to extract start and end positions."""
    # Remove any non-digit, non-hyphen characters
    clean_region = re.sub(r'[^\d-]', '', region_str)
    
    if '-' in clean_region:
        parts = clean_region.split('-')
        if len(parts) == 2:
            try:
                start = int(parts[0])
                end = int(parts[1])
                return (start, end)
            except ValueError:
                pass
    
    return None

def correct_primer_positions(df, consensus_seqs, min_match=0.8, min_diff=5):
    """Check and correct primer binding positions."""
    # Make a copy of the dataframe to modify
    corrected_df = df.copy()
    
    # Track which rows were corrected
    corrections = []
    
    # Process each row
    for idx, row in df.iterrows():
        # Determine which segment this primer is for
        if row['F Name'].startswith('S'):
            segment = 'S'
        elif row['F Name'].startswith('M'):
            segment = 'M'
        else:
            # Skip primers for other segments
            continue
        
        consensus_seq = consensus_seqs[segment]
        
        # Extract primer sequences
        f_seq_col = [col for col in df.columns if '5\'-3\' Sequence' in col and col.startswith('F')][0]
        r_seq_col = [col for col in df.columns if '5\'-3\' Sequence' in col and col.startswith('R')][0]
        
        f_seq = row[f_seq_col]
        r_seq = row[r_seq_col]
        
        # Extract current positions
        f_region = row.get('F Region', '')
        r_region = row.get('R Region', '')
        
        # Parse current positions
        current_f_pos = parse_region(f_region)
        current_r_pos = parse_region(r_region)
        
        # Check forward primer binding
        new_f_pos = verify_primer_binding(f_seq, consensus_seq, is_reverse=False, min_match=min_match)
        
        # Check reverse primer binding
        new_r_pos = verify_primer_binding(r_seq, consensus_seq, is_reverse=True, min_match=min_match)
        
        # Initialize correction info with accurate original positions
        correction_info = {
            'F Name': row['F Name'],
            'R Name': row['Reverse Name'],
            'F Old Position': f_region,
            'R Old Position': r_region,
            'F Corrected': False,
            'R Corrected': False,
            'F Match Score': new_f_pos[2] if new_f_pos else 0,
            'R Match Score': new_r_pos[2] if new_r_pos else 0,
            'F New Position': f_region,  # Default to original
            'R New Position': r_region   # Default to original
        }
        
        # Only update if we have high confidence and position is significantly different
        if new_f_pos and current_f_pos and (
                abs(new_f_pos[0] - current_f_pos[0]) > min_diff or 
                abs(new_f_pos[1] - current_f_pos[1]) > min_diff):
            # Update the F Region
            new_f_region = update_region_format(f_region, new_f_pos[0], new_f_pos[1], is_reverse=False)
            corrected_df.at[idx, 'F Region'] = new_f_region
            correction_info['F New Position'] = new_f_region
            correction_info['F Corrected'] = True
        
        if new_r_pos and current_r_pos and (
                abs(new_r_pos[0] - current_r_pos[0]) > min_diff or 
                abs(new_r_pos[1] - current_r_pos[1]) > min_diff):
            # Update the R Region
            new_r_region = update_region_format(r_region, new_r_pos[0], new_r_pos[1], is_reverse=True)
            corrected_df.at[idx, 'R Region'] = new_r_region
            correction_info['R New Position'] = new_r_region
            correction_info['R Corrected'] = True
        
        # Update Overall Region if needed - use amplicon boundaries
        if (correction_info['F Corrected'] or correction_info['R Corrected']) and 'Overall Region' in df.columns:
            # Calculate new overall region using actual amplicon boundaries
            f_start = parse_region(correction_info['F New Position'])[0] if parse_region(correction_info['F New Position']) else None
            r_end = parse_region(correction_info['R New Position'])[1] if parse_region(correction_info['R New Position']) else None
            
            if f_start is not None and r_end is not None:
                corrected_df.at[idx, 'Overall Region'] = f"{f_start}-{r_end}"
                correction_info['Overall Region Updated'] = True
        
        # Add to corrections list if changed
        if correction_info['F Corrected'] or correction_info['R Corrected']:
            corrections.append(correction_info)
    
    return corrected_df, corrections

def main():
    args = parse_args()
    
    # Load primers
    df = load_primers(args.primers)
    
    # Load consensus sequences
    consensus_seqs = load_consensus_sequences(args.s_consensus, args.m_consensus)
    
    # Correct primer positions
    print("Checking primer binding positions...")
    corrected_df, corrections = correct_primer_positions(
        df, consensus_seqs, min_match=args.min_match, min_diff=args.min_diff
    )
    
    # Save corrected CSV
    corrected_df.to_csv(args.output, index=False, encoding='utf-8')
    print(f"Corrected primers saved to {args.output}")
    
    # Save corrections summary
    if corrections:
        correction_summary_path = args.output.replace('.csv', '_corrections.csv')
        corrections_df = pd.DataFrame(corrections)
        corrections_df.to_csv(correction_summary_path, index=False, encoding='utf-8')
        print(f"Corrections summary saved to {correction_summary_path}")
        
        # Print summary
        print(f"\nCorrected {len(corrections)} primer pairs:")
        for c in corrections:
            f_change = f"{c['F Old Position']} → {c['F New Position']}" if c['F Corrected'] else "No change"
            r_change = f"{c['R Old Position']} → {c['R New Position']}" if c['R Corrected'] else "No change"
            print(f"  - {c['F Name']}/{c['R Name']}: Forward: {f_change}, Reverse: {r_change}")
    else:
        print("No corrections needed.")

if __name__ == "__main__":
    main() 