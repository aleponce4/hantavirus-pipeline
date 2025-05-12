#!/usr/bin/env python3
"""
Primer Evaluation Script for Hantavirus Pipeline

This script evaluates primers for quality based on various criteria:
- GC Content
- Melting Temperature
- Secondary Structure (Hairpins)
- BLAST E-value
- Sequence length and composition
"""

import pandas as pd
import os
import sys
import argparse
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio import SeqIO
import subprocess
import seqfold
from seqfold import fold, dg, dg_cache, dot_bracket
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
from dna_features_viewer import GraphicFeature, GraphicRecord
import numpy as np
import glob
import pysam  # For BAM file processing

# Base directories
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS_DIR = os.path.join(BASE_DIR, 'results', 'primer_evaluation')
os.makedirs(RESULTS_DIR, exist_ok=True)

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Evaluate primers for quality.')
    parser.add_argument('--primers', required=True, help='Path to primers CSV file')
    parser.add_argument('--segment', required=True, choices=['S_segment', 'M_segment', 'L_segment'], 
                        help='Viral segment to analyze')
    parser.add_argument('--blast_db', required=True, help='Path to BLAST database')
    parser.add_argument('--consensus', required=True, help='Path to consensus sequence file')
    parser.add_argument('--output', help='Output directory', default=RESULTS_DIR)
    
    return parser.parse_args()

def load_primers(primers_file):
    """Load primers from CSV file."""
    print(f"Loading primers from {primers_file}")
    try:
        # Use UTF-8 encoding when reading the CSV file
        df = pd.read_csv(primers_file, encoding='utf-8')
        
        # Clean column names by stripping whitespace and replacing special characters
        df.columns = df.columns.str.strip()
        clean_col_names = {}
        for col in df.columns:
            clean_name = col.replace('´', "'").replace('°C', ' C').replace('ºC', ' C')
            if clean_name != col:
                clean_col_names[col] = clean_name
        
        # Rename columns if needed
        if clean_col_names:
            df = df.rename(columns=clean_col_names)
        
        # Clean up text fields to remove non-ASCII characters
        for col in df.columns:
            if df[col].dtype == 'object':  # Only process string columns
                # Replace special characters like the degree symbol (°) with plain text
                df[col] = df[col].str.replace('°C', ' C', regex=False) if isinstance(df[col], pd.Series) else df[col]
                df[col] = df[col].str.replace('ºC', ' C', regex=False) if isinstance(df[col], pd.Series) else df[col]
                # Replace other problematic characters
                df[col] = df[col].str.replace('´', "'", regex=False) if isinstance(df[col], pd.Series) else df[col]
                
        # Filter only primers for the specified segment
        return df
    except Exception as e:
        print(f"Error loading primers: {e}")
        sys.exit(1)

def calculate_primer_properties(primer_sequence):
    """Calculate GC content and melting temperature for a primer."""
    # Clean the sequence (remove spaces or any non-nucleotide characters)
    cleaned_seq = ''.join(c for c in primer_sequence if c.upper() in 'ACGT')
    
    if not cleaned_seq:
        return 0, 0  # Return zeros for empty or invalid sequences
    
    # Calculate GC content using gc_fraction instead of GC
    gc_content = gc_fraction(cleaned_seq)
    
    # Calculate melting temperature using Biopython's MeltingTemp module
    # Updated parameters to match standard PCR conditions:
    # - 1.5 mM Mg²⁺
    # - 50 mM Na⁺/K⁺
    # - 0.6 mM dNTPs
    melting_temp = mt.Tm_NN(
        Seq(cleaned_seq),
        dnac1=250,  # nanomolar
        dnac2=250,  # nanomolar
        Na=0,      # 50 mM Na⁺
        K=50,       # 50 mM K⁺
        Tris=0,
        Mg=1.5,     # 1.5 mM Mg²⁺
        dNTPs=0.6,  # 0.6 mM dNTPs
        saltcorr=7,
        nn_table=mt.DNA_NN3,  # Thermodynamic NN values from Allawi and SantaLucia 1997
    )
    
    return gc_content, melting_temp

def predict_gibbs_free_energy(primer_sequence):
    """Predict Gibbs Free Energy of the most stable hairpin."""
    # Clean the sequence
    cleaned_seq = ''.join(c for c in primer_sequence if c.upper() in 'ACGT')
    
    if not cleaned_seq:
        return 0  # Return zero for empty or invalid sequences
    
    # Predict Gibbs Free Energy
    try:
        gibb_free_energy = dg(cleaned_seq, temp=20.0)
        # Replace positive values more than 10 with 2 (as per the original code)
        if gibb_free_energy > 10:
            gibb_free_energy = 2
        return gibb_free_energy
    except Exception as e:
        print(f"Error calculating Gibbs free energy: {e}")
        return 0

def get_reverse_complement(seq):
    """Generate the reverse complement of a sequence."""
    return str(Seq(seq).reverse_complement())

def run_blast_and_get_evalue(primer_sequence, db_path):
    """Run BLAST and get the E-value for a primer."""
    # Clean the sequence
    cleaned_seq = ''.join(c for c in primer_sequence if c.upper() in 'ACGT')
    
    if not cleaned_seq:
        return None  # Return None for empty or invalid sequences
    
    # Write the sequence to a temporary file directly (no reverse complement)
    temp_dir = os.path.join(RESULTS_DIR, 'temp')
    os.makedirs(temp_dir, exist_ok=True)
    temp_file = os.path.join(temp_dir, "temp_primer.fasta")
    
    # Write the primer sequence to a temporary file
    with open(temp_file, "w") as f:
        f.write(f">primer\n{cleaned_seq}\n")
    
    try:
        # Use more sensitive BLAST parameters for short primers
        cmd = [
            "blastn", "-query", temp_file, "-db", db_path,
            "-outfmt", "6",
            "-evalue", "10",         # Allow higher E-value threshold
            "-word_size", "7",        # Smaller word size for increased sensitivity with short sequences
            "-task", "blastn-short",  # Use blastn-short for short sequences
            "-dust", "no"             # Disable dust filtering
        ]
        
        print(f"Running BLAST for: {cleaned_seq}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Check if BLAST ran successfully
        if result.returncode != 0:
            print(f"Error running BLAST: {result.stderr}")
            return None

        # Parse the BLAST output to get the E-value
        lines = result.stdout.strip().split('\n')
        if lines and lines[0]:
            try:
                # Extract the E-value from the first hit
                evalue = float(lines[0].split('\t')[-2])
                return evalue
            except (IndexError, ValueError) as e:
                print(f"Error parsing BLAST output: {e}")
                return None
        else:
            # No hits found, helpful debug output
            print(f"No BLAST hits found for: {cleaned_seq}")
            return None
    except Exception as e:
        print(f"Error running BLAST: {e}")
        return None
    finally:
        # Clean up temporary file
        if os.path.exists(temp_file):
            os.remove(temp_file)

def verify_primer_binding_site(primer_seq, reference_seq, position, is_reverse=False, flanking_region=5):
    """
    Verify if a primer matches its intended binding region in the reference sequence.
    
    Args:
        primer_seq (str): Primer sequence
        reference_seq (str): Full reference sequence
        position (tuple): Tuple containing (start, end) positions of primer
        is_reverse (bool): Whether this is a reverse primer
        flanking_region (int): Number of bases to check on either side of the expected position
    
    Returns:
        dict: Dictionary with match percentage, actual position, and match status
    """
    # Handle None position by searching the entire sequence
    if position is None:
        print(f"Warning: No position provided for {'reverse' if is_reverse else 'forward'} primer. Searching entire sequence.")
        # Slide the primer across the entire reference sequence
        primer_to_match = primer_seq
        if is_reverse:
            primer_to_match = str(Seq(primer_seq).reverse_complement())
        
        best_match = 0
        best_pos = None
        
        # Scan the entire sequence for the best match
        for i in range(len(reference_seq) - len(primer_to_match) + 1):
            match_count = sum(a == b for a, b in zip(primer_to_match, reference_seq[i:i+len(primer_to_match)]))
            match_percent = match_count / len(primer_to_match) * 100
            
            if match_percent > best_match:
                best_match = match_percent
                best_pos = i
        
        match_status = "Good" if best_match >= 80 else "Poor"
        return {
            "match_percent": best_match,
            "actual_position": (best_pos, best_pos + len(primer_to_match) - 1) if best_pos is not None else (None, None),
            "status": match_status,
            "distance_from_expected": None  # No expected position
        }
    
    # Normal case where position is provided
    start, end = position
    
    # Adjust for 0-based indexing
    start = max(0, start)
    end = min(len(reference_seq) - 1, end)
    
    # Expand region to check for slight positional variations
    check_start = max(0, start - flanking_region)
    check_end = min(len(reference_seq), end + flanking_region + len(primer_seq))
    
    # Extract the region to check
    region_to_check = reference_seq[check_start:check_end]
    
    # For reverse primers, we need to use the reverse complement
    primer_to_match = primer_seq
    if is_reverse:
        primer_to_match = str(Seq(primer_seq).reverse_complement())
    
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
    
    return {
        "match_percent": best_match,
        "actual_position": (best_pos, best_pos + len(primer_to_match) - 1) if best_pos is not None else (None, None),
        "status": match_status,
        "distance_from_expected": best_pos - start if best_pos is not None else None
    }

def check_primer_quality(row):
    """Check primer quality against established criteria."""
    problems = []
    
    # Extract primer sequences
    f_seq = row.get('F 5\'-3\' Sequence', '') if 'F 5\'-3\' Sequence' in row else row.get('F 5´-3´ Sequence', '')
    r_seq = row.get('R 5\'-3\' Sequence', '') if 'R 5\'-3\' Sequence' in row else row.get('R 5´-3´ Sequence', '')
    
    # Skip if primers are missing
    if not f_seq or not r_seq:
        return "Missing primer sequence"
    
    # Clean the sequences
    f_seq_clean = ''.join(c for c in f_seq if c.upper() in 'ACGT')
    r_seq_clean = ''.join(c for c in r_seq if c.upper() in 'ACGT')
    
    # Avoid long runs of the same bases or repeats of consecutive dinucleotides
    repeats = ['AAAAA', 'TTTTT', 'CCCCC', 'GGGGG', 'ATATAT', 'TATATA', 'CGCGCG', 'GCGCGC']
    if any(rep in f_seq_clean for rep in repeats):
        problems.append('Long run of bases or dinucleotide repeats in Forward Primer')
    if any(rep in r_seq_clean for rep in repeats):
        problems.append('Long run of bases or dinucleotide repeats in Reverse Primer')

    # Check G/C content
    if not (0.35 <= row['F_GC_Content'] <= 0.65):
        problems.append(f'Forward Primer G/C content out of range: {row["F_GC_Content"]:.2f}')
    if not (0.35 <= row['R_GC_Content'] <= 0.65):
        problems.append(f'Reverse Primer G/C content out of range: {row["R_GC_Content"]:.2f}')
  
    # Check melting temperature
    if not (58 <= row['F_Melting_Temperature'] <= 72):
        problems.append(f'Forward Primer Tm out of range: {row["F_Melting_Temperature"]:.2f} C')
    if not (58 <= row['R_Melting_Temperature'] <= 72):
        problems.append(f'Reverse Primer Tm out of range: {row["R_Melting_Temperature"]:.2f} C')
    
    # Check primer length
    if not (18 <= len(f_seq_clean) <= 35):
        problems.append(f'Forward Primer length out of range: {len(f_seq_clean)} bp')
    if not (18 <= len(r_seq_clean) <= 35):
        problems.append(f'Reverse Primer length out of range: {len(r_seq_clean)} bp')
    
    # Check E-value
    if pd.notna(row['F_BLAST_E_Value']) and row['F_BLAST_E_Value'] > 1e-4:
        problems.append(f'Forward Primer BLAST E-value too high: {row["F_BLAST_E_Value"]}')
    if pd.notna(row['R_BLAST_E_Value']) and row['R_BLAST_E_Value'] > 1e-4:
        problems.append(f'Reverse Primer BLAST E-value too high: {row["R_BLAST_E_Value"]}')
    
    # Check Gibbs Free Energy
    if row['F_Gibbs_Free_Energy'] < -6:
        problems.append(f'Forward Primer Gibbs Free Energy too low: {row["F_Gibbs_Free_Energy"]:.2f} kcal/mol')
    if row['R_Gibbs_Free_Energy'] < -6:
        problems.append(f'Reverse Primer Gibbs Free Energy too low: {row["R_Gibbs_Free_Energy"]:.2f} kcal/mol')
    
    # Check binding site match if data is available
    if 'F_Binding_Match_Percent' in row and row['F_Binding_Match_Percent'] < 80:
        problems.append(f'Forward Primer binding site match too low: {row["F_Binding_Match_Percent"]:.1f}%')
    if 'R_Binding_Match_Percent' in row and row['R_Binding_Match_Percent'] < 80:
        problems.append(f'Reverse Primer binding site match too low: {row["R_Binding_Match_Percent"]:.1f}%')
    
    return ', '.join(problems) if problems else 'Pass'

def get_coverage_from_bam(bam_file, min_quality=0, max_coverage=1000):
    """Extract coverage data from a BAM file"""
    coverage = []
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            # Get reference length from the first reference
            ref_name = bam.references[0]
            ref_length = bam.lengths[0]
            
            # Generate coverage array
            coverage = np.zeros(ref_length)
            for pileupcol in bam.pileup(ref_name, min_mapping_quality=min_quality):
                # Cap coverage at max_coverage value
                coverage[pileupcol.pos] = min(pileupcol.n, max_coverage)
    except Exception as e:
        print(f"Error processing BAM file {bam_file}: {e}")
    
    return coverage

def find_bam_files(segment):
    """Find all BAM files for a specific segment"""
    bam_files = []
    
    # Look in results directory for sample subdirectories
    results_path = os.path.join(BASE_DIR, 'results')
    for sample_dir in os.listdir(results_path):
        # Skip special directories
        if sample_dir in ['plots', 'first_pass', 'primer_evaluation'] or not os.path.isdir(os.path.join(results_path, sample_dir)):
            continue
        
        # Check for alignment directory for this segment
        align_dir = os.path.join(results_path, sample_dir, 'alignment', segment)
        if not os.path.isdir(align_dir):
            continue
        
        # Look for BAM file
        bam_file = os.path.join(align_dir, f"{sample_dir}.bam")
        if os.path.exists(bam_file):
            bam_files.append(bam_file)
    
    return bam_files

def get_average_coverage(segment, min_quality=0, max_coverage=1000):
    """Calculate average coverage across all samples for a segment"""
    bam_files = find_bam_files(segment)
    
    if not bam_files:
        print(f"No BAM files found for {segment}")
        return np.array([])
    
    # Get coverage from each BAM file
    all_coverages = []
    for bam_file in bam_files:
        coverage = get_coverage_from_bam(bam_file, min_quality, max_coverage)
        if len(coverage) > 0:
            all_coverages.append(coverage)
    
    if not all_coverages:
        return np.array([])
    
    # Pad shorter arrays with zeros to match the longest
    max_length = max(len(cov) for cov in all_coverages)
    padded_coverages = []
    for cov in all_coverages:
        if len(cov) < max_length:
            padded = np.zeros(max_length)
            padded[:len(cov)] = cov
            padded_coverages.append(padded)
        else:
            padded_coverages.append(cov)
    
    # Calculate average coverage
    average_coverage = np.mean(padded_coverages, axis=0)
    
    return average_coverage

def calculate_amplicon_coverage(f_start, r_end, segment):
    """Calculate the average coverage for an amplicon region"""
    # Get average coverage across all samples
    average_coverage = get_average_coverage(segment)
    
    if len(average_coverage) == 0:
        return 0.0
    
    # Ensure indices are within valid bounds
    f_start = max(0, min(f_start, len(average_coverage)-1))
    r_end = max(0, min(r_end, len(average_coverage)-1))
    
    # Extract coverage for the amplicon region
    amplicon_coverage = average_coverage[f_start:r_end+1]
    
    # Calculate average
    if len(amplicon_coverage) > 0:
        return np.mean(amplicon_coverage)
    else:
        return 0.0

def calculate_amplicon_dropout(f_start, r_end, segment, threshold=10):
    """Calculate the percentage of bases within an amplicon that fall below the coverage threshold."""
    # Get average coverage across all samples
    average_coverage = get_average_coverage(segment)
    
    if len(average_coverage) == 0:
        return 100.0  # Return 100% dropout if no coverage data
    
    # Ensure indices are within valid bounds
    f_start = max(0, min(f_start, len(average_coverage)-1))
    r_end = max(0, min(r_end, len(average_coverage)-1))
    
    # Extract coverage for the amplicon region
    amplicon_coverage = average_coverage[f_start:r_end+1]
    
    if len(amplicon_coverage) == 0:
        return 100.0
    
    # Count bases below threshold
    bases_below_threshold = np.sum(amplicon_coverage < threshold)
    
    # Calculate dropout percentage
    dropout_percentage = (bases_below_threshold / len(amplicon_coverage)) * 100
    
    return dropout_percentage

def plot_primers_to_target(df, consensus_seq, output_path, segment):
    """Plot primers relative to the target consensus sequence with coverage plot"""
    # Get the consensus sequence length
    consensus_sequence_length = len(consensus_seq)
    
    # For M segment, we need to extend the plot to cover the UTR region up to 4205
    # The CDS ends at position 3682
    if segment.startswith('M'):
        plot_sequence_length = 4205
        cds_end = 3682
    else:
        plot_sequence_length = consensus_sequence_length
        cds_end = None
    
    # Adjust figure width based on segment size
    # Typical segment lengths: S ~1.8kb, M ~3.6kb, L ~6.5kb
    # Calculate a segment-dependent width multiplier
    if segment.startswith('S'):
        width_multiplier = 1.0
    elif segment.startswith('M'):
        width_multiplier = 1.5
    elif segment.startswith('L'):
        width_multiplier = 2.0
    else:
        width_multiplier = 1.0
    
    # Base width starts at 15 inches
    base_width = 15
    adjusted_width = base_width * width_multiplier
    
    # Get average coverage data
    average_coverage = get_average_coverage(segment)
    
    # If no coverage data available, just do the regular primer plot
    if len(average_coverage) == 0:
        has_coverage = False
        fig_height = 3 + 0.2 * len(df)  # Reduced height
        fig, ax = plt.subplots(1, figsize=(adjusted_width, fig_height))  # Segment-dependent width
        
        # Create features for the original sequence - using the standard approach
        # but with better organization
        features = [
            GraphicFeature(start=0, end=plot_sequence_length, strand=0, color="#ffffff", label="Consensus")
        ]
        
        # Create a color map for different primer pairs
        primer_colors = plt.cm.tab10.colors
        
        # First add all amplicon features
        for i, row in df.iterrows():
            color = primer_colors[i % len(primer_colors)]
            
            try:
                # Extract positions from the region fields
                f_region = row.get('F Region', '')
                r_region = row.get('R Region', '')
                
                # Parse position from region field (e.g., "1-32F" -> start=1, end=32)
                f_match = f_region.split('-')
                f_start = int(f_match[0])
                f_end = int(''.join(c for c in f_match[1] if c.isdigit()))
                
                r_match = r_region.split('-')
                r_start = int(r_match[0])
                r_end = int(''.join(c for c in r_match[1] if c.isdigit()))
                
                # Adjust 1-based positions to 0-based
                f_start -= 1
                f_end -= 1
                r_start -= 1
                r_end -= 1
                
                # Calculate amplicon length
                amplicon_length = r_end - f_start + 1
                
                # Create a more compact label format
                if amplicon_length >= 1000:
                    size_label = f"{amplicon_length/1000:.1f}kb"
                else:
                    size_label = f"{amplicon_length}bp"
                
                # Create amplicon feature (will be displayed in top lanes)
                features.append(
                    GraphicFeature(start=f_start, end=r_end, strand=0, color=color, thickness=12, 
                                  label=f"A{i+1}({size_label})", fontsize=5, label_position="middle", 
                                  box_color=color, text_color="white")
                )
            except Exception as e:
                print(f"Error adding amplicon feature for pair {i+1}: {e}")
        
        # Then add all primer features
        for i, row in df.iterrows():
            color = primer_colors[i % len(primer_colors)]
            
            try:
                f_name = row.get('F Name', '')
                r_name = row.get('Reverse Name', '')
                
                # Extract positions from the region fields
                f_region = row.get('F Region', '')
                r_region = row.get('R Region', '')
                
                # Parse position from region field
                f_match = f_region.split('-')
                f_start = int(f_match[0])
                f_end = int(''.join(c for c in f_match[1] if c.isdigit()))
                
                r_match = r_region.split('-')
                r_start = int(r_match[0])
                r_end = int(''.join(c for c in r_match[1] if c.isdigit()))
                
                # Adjust 1-based positions to 0-based
                f_start -= 1
                f_end -= 1
                r_start -= 1
                r_end -= 1
                
                # Add forward primer and reverse primer features
                features.append(
                    GraphicFeature(start=f_start, end=f_end, strand=1, color=color, label=f_name, fontsize=5)
                )
                
                features.append(
                    GraphicFeature(start=r_start, end=r_end, strand=-1, color=color, label=r_name, fontsize=5)
                )
            except Exception as e:
                print(f"Error adding primer features for pair {i+1}: {e}")
        
        # Create the GraphicRecord
        record = GraphicRecord(sequence_length=plot_sequence_length, features=features)
        
        # Plot the GraphicRecord
        record.plot(ax=ax)
        ax.set_title(f"Primer Pairs and Amplicons on {segment} Consensus Sequence")
        
        # Mark the CDS end with a vertical line for M segment
        if cds_end is not None:
            ax.axvline(x=cds_end, color='red', linestyle='--', linewidth=1)
            ax.text(cds_end + 10, 0.95, "3' UTR →", transform=ax.get_xaxis_transform(), 
                    color='red', fontsize=8, fontweight='bold')
        
    else:
        has_coverage = True
        # To maintain consistency, always get max_coverage
        max_coverage = max(average_coverage) * 1.1 if average_coverage else 100
        
        # Setup a more complex figure with 2 subplots for coverage and primers
        fig_height = 5  # Fixed height for the combined plot
        
        # Create figure and gridspec for more control over subplot sizes
        fig = plt.figure(figsize=(adjusted_width, fig_height))
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 2], hspace=0.1)
        
        # Set up coverage plot in the top panel
        ax_coverage = fig.add_subplot(gs[0])
        
        # Plot coverage
        ax_coverage.fill_between(
            range(1, len(average_coverage) + 1),
            average_coverage,
            color="lightgrey",
            alpha=0.8,
            label="Coverage"
        )
        
        # Add horizontal line at the low coverage threshold
        LOW_COVERAGE = 10
        ax_coverage.axhline(y=LOW_COVERAGE, color='red', linestyle='--', 
                           linewidth=0.8, alpha=0.7, label=f"Low coverage (<{LOW_COVERAGE}x)")
        
        # Highlight low coverage regions more prominently
        low_coverage_regions = np.array(average_coverage) < LOW_COVERAGE
        if np.any(low_coverage_regions):
            ax_coverage.fill_between(
                range(1, len(average_coverage) + 1),
                [LOW_COVERAGE if x < LOW_COVERAGE else 0 for x in average_coverage],
                color="red",
                alpha=0.3,
                label="Dropout regions"
            )
        
        # Add legend
        ax_coverage.legend(loc='upper right', fontsize=6)
        
        # Format the coverage plot
        ax_coverage.set_title(f"Coverage and Primers - {segment}", fontsize=12)
        ax_coverage.set_ylabel("Coverage", fontsize=8)
        
        # Clean up the coverage plot
        for spine in ['top', 'right']:
            ax_coverage.spines[spine].set_visible(False)
            
        # Set y-axis limit to keep low coverage visible but also see high coverage peaks
        ax_coverage.set_ylim([0, max_coverage])
        ax_coverage.set_xlim([1, plot_sequence_length])
        ax_coverage.set_xlabel('')
        ax_coverage.set_xlim([1, plot_sequence_length])
        ax_primers.set_xlim([1, plot_sequence_length])
        
        # Mark the CDS end with a vertical line for M segment
        if cds_end is not None:
            ax_coverage.axvline(x=cds_end, color='red', linestyle='--', linewidth=1)
            ax_coverage.text(cds_end + 10, max_coverage * 0.9, "3' UTR →", 
                    color='red', fontsize=8, fontweight='bold')
            ax_primers.axvline(x=cds_end, color='red', linestyle='--', linewidth=1)
            ax_primers.text(cds_end + 10, 0.95, "3' UTR →", transform=ax_primers.get_xaxis_transform(), 
                    color='red', fontsize=8, fontweight='bold')
        
        # Clean up spines
        for spine in ['top', 'right']:
            ax_primers.spines[spine].set_visible(False)
    
    # Final formatting and save
    plt.tight_layout(h_pad=0.5)  # Reduce padding between subplots
    plt.savefig(output_path, dpi=300, bbox_inches='tight')  # Crop any extra whitespace
    
    # Also save as SVG for vector editing
    svg_output_path = output_path.replace('.png', '.svg')
    plt.savefig(svg_output_path, format='svg', bbox_inches='tight')
    plt.close()
    
    # Now create a second plot with JUST the amplicons (no primer arrows and labels)
    plot_amplicons_only(df, consensus_seq, output_path.replace('.png', '_amplicons_only.png'), segment, plot_sequence_length, cds_end)

def plot_amplicons_only(df, consensus_seq, output_path, segment, plot_sequence_length, cds_end=None):
    """Create a simplified plot with just the amplicons, no primer arrows or labels"""
    
    # Adjust figure width based on segment size
    if segment.startswith('S'):
        width_multiplier = 1.0
    elif segment.startswith('M'):
        width_multiplier = 1.5
    elif segment.startswith('L'):
        width_multiplier = 2.0
    else:
        width_multiplier = 1.0
    
    # Base width starts at 15 inches
    base_width = 15
    adjusted_width = base_width * width_multiplier
    
    # Setup the figure
    fig_height = 3  # Reduced height for amplicons only
    fig, ax = plt.subplots(1, figsize=(adjusted_width, fig_height))
    
    # Create features for the consensus sequence
    features = [
        GraphicFeature(start=0, end=plot_sequence_length, strand=0, color="#ffffff", label="Consensus")
    ]
    
    # Create a color map for different primer pairs
    primer_colors = plt.cm.tab10.colors
    
    # Add all amplicon features
    for i, row in df.iterrows():
        color = primer_colors[i % len(primer_colors)]
        
        try:
            # Extract positions from the region fields
            f_region = row.get('F Region', '')
            r_region = row.get('R Region', '')
            
            # Parse position from region field
            f_match = f_region.split('-')
            f_start = int(f_match[0])
            f_end = int(''.join(c for c in f_match[1] if c.isdigit()))
            
            r_match = r_region.split('-')
            r_start = int(r_match[0])
            r_end = int(''.join(c for c in r_match[1] if c.isdigit()))
            
            # Adjust 1-based positions to 0-based
            f_start -= 1
            f_end -= 1
            r_start -= 1
            r_end -= 1
            
            # Calculate amplicon length
            amplicon_length = r_end - f_start + 1
            
            # Create a more compact label format
            if amplicon_length >= 1000:
                size_label = f"{amplicon_length/1000:.1f}kb"
            else:
                size_label = f"{amplicon_length}bp"
            
            # Create amplicon feature with increased thickness for better visibility
            features.append(
                GraphicFeature(start=f_start, end=r_end, strand=0, color=color, thickness=14, 
                              label=f"A{i+1}({size_label})", fontsize=5, label_position="middle", 
                              box_color=color, text_color="white")
            )
        except Exception as e:
            print(f"Error adding amplicon feature for pair {i+1}: {e}")
    
    # Create and plot the GraphicRecord
    record = GraphicRecord(sequence_length=plot_sequence_length, features=features)
    record.plot(ax=ax)
    
    # Set title and format
    ax.set_title(f"Amplicons Only - {segment} Consensus Sequence")
    
    # Mark the CDS end with a vertical line for M segment
    if cds_end is not None:
        ax.axvline(x=cds_end, color='red', linestyle='--', linewidth=1)
        ax.text(cds_end + 10, 0.95, "3' UTR →", transform=ax.get_xaxis_transform(), 
                color='red', fontsize=8, fontweight='bold')
    
    # Clean up spines
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    
    # Save the plot
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    # Also save as SVG
    svg_output_path = output_path.replace('.png', '.svg')
    plt.savefig(svg_output_path, format='svg', bbox_inches='tight')
    
    plt.close()

def extract_positions(row):
    """Extract forward and reverse primer positions to determine amplicon boundaries"""
    try:
        # Extract positions from the region fields
        f_region = row.get('F Region', '')
        r_region = row.get('R Region', '')
        
        if not f_region or not r_region:
            return None
        
        # Parse positions from region field (e.g., "1-32F" -> start=1, end=32)
        # Handle potential formatting issues
        try:
            f_match = f_region.split('-')
            if len(f_match) != 2:
                print(f"Warning: Invalid forward region format: {f_region}")
                return None
                
            f_start = int(f_match[0])
            f_end = int(''.join(c for c in f_match[1] if c.isdigit()))
        except (ValueError, IndexError) as e:
            print(f"Warning: Could not parse forward region {f_region}: {e}")
            return None
            
        try:
            r_match = r_region.split('-')
            if len(r_match) != 2:
                print(f"Warning: Invalid reverse region format: {r_region}")
                return None
                
            r_start = int(r_match[0])
            r_end = int(''.join(c for c in r_match[1] if c.isdigit()))
        except (ValueError, IndexError) as e:
            print(f"Warning: Could not parse reverse region {r_region}: {e}")
            return None
        
        # Adjust to 0-based positions for internal calculations
        f_start -= 1
        r_end -= 1
        
        return (f_start, r_end)
    except Exception as e:
        print(f"Error extracting positions: {e}")
        return None

def main():
    args = parse_args()
    
    # Load primers
    df = load_primers(args.primers)
    
    # Filter primers for the specified segment
    segment_prefix = args.segment[0]  # S, M, or L
    df = df[df['F Name'].str.startswith(segment_prefix)]
    
    if df.empty:
        print(f"No primers found for segment {args.segment}")
        sys.exit(1)
    
    print(f"Loaded {len(df)} primer pairs for {args.segment}")
    
    # Determine the sequence column names
    f_seq_col = [col for col in df.columns if '5\'-3\' Sequence' in col and col.startswith('F')][0]
    r_seq_col = [col for col in df.columns if '5\'-3\' Sequence' in col and col.startswith('R')][0]
    
    # Load reference sequence
    consensus_seq = str(SeqIO.read(args.consensus, "fasta").seq)
    
    # Calculate properties for each primer
    print("Calculating primer properties...")
    
    # Forward primers
    df['F_GC_Content'] = df[f_seq_col].apply(lambda x: calculate_primer_properties(x)[0])
    df['F_Melting_Temperature'] = df[f_seq_col].apply(lambda x: calculate_primer_properties(x)[1])
    df['F_Gibbs_Free_Energy'] = df[f_seq_col].apply(predict_gibbs_free_energy)
    
    # Reverse primers  
    df['R_GC_Content'] = df[r_seq_col].apply(lambda x: calculate_primer_properties(x)[0])
    df['R_Melting_Temperature'] = df[r_seq_col].apply(lambda x: calculate_primer_properties(x)[1])
    df['R_Gibbs_Free_Energy'] = df[r_seq_col].apply(predict_gibbs_free_energy)
    
    # Run BLAST for each primer
    print("Running BLAST analysis...")
    df['F_BLAST_E_Value'] = df[f_seq_col].apply(lambda x: run_blast_and_get_evalue(x, args.blast_db))
    df['R_BLAST_E_Value'] = df[r_seq_col].apply(lambda x: run_blast_and_get_evalue(x, args.blast_db))
    
    # Calculate average coverage and dropout for each amplicon
    print("Calculating amplicon coverage and dropout...")
    df['Amplicon_Positions'] = df.apply(lambda row: extract_positions(row), axis=1)
    
    # Verify primer binding site match
    print("Verifying primer binding sites...")
    
    # Process each row to get the expected primer positions
    for idx, row in df.iterrows():
        f_region = row.get('F Region', '')
        r_region = row.get('R Region', '')
        
        # Parse positions
        try:
            # Setup default values
            f_start = None
            f_end = None
            r_start = None
            r_end = None
            
            # Forward primer
            try:
                if f_region and '-' in f_region:
                    f_match = f_region.split('-')
                    if len(f_match) == 2:
                        f_start = int(f_match[0]) - 1  # Convert to 0-based
                        f_end = int(''.join(c for c in f_match[1] if c.isdigit())) - 1  # Convert to 0-based
            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse forward region {f_region} for row {idx}: {e}")
            
            # Reverse primer
            try:
                if r_region and '-' in r_region:
                    r_match = r_region.split('-')
                    if len(r_match) == 2:
                        r_start = int(r_match[0]) - 1  # Convert to 0-based
                        r_end = int(''.join(c for c in r_match[1] if c.isdigit())) - 1  # Convert to 0-based
            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse reverse region {r_region} for row {idx}: {e}")
            
            # Verify binding - handle None positions by searching entire sequence
            f_binding = verify_primer_binding_site(
                row[f_seq_col], 
                consensus_seq, 
                (f_start, f_end) if f_start is not None and f_end is not None else None, 
                is_reverse=False
            )
            
            r_binding = verify_primer_binding_site(
                row[r_seq_col], 
                consensus_seq, 
                (r_start, r_end) if r_start is not None and r_end is not None else None, 
                is_reverse=True
            )
            
            # Store results
            df.at[idx, 'F_Binding_Match_Percent'] = f_binding['match_percent']
            df.at[idx, 'F_Binding_Status'] = f_binding['status']
            df.at[idx, 'F_Binding_Offset'] = f_binding['distance_from_expected']
            
            df.at[idx, 'R_Binding_Match_Percent'] = r_binding['match_percent']
            df.at[idx, 'R_Binding_Status'] = r_binding['status']
            df.at[idx, 'R_Binding_Offset'] = r_binding['distance_from_expected']
            
        except Exception as e:
            print(f"Error verifying primer binding for row {idx}: {e}")
            df.at[idx, 'F_Binding_Status'] = "Error"
            df.at[idx, 'R_Binding_Status'] = "Error"
    
    df['Average_Amplicon_Coverage'] = df.apply(
        lambda row: calculate_amplicon_coverage(
            row['Amplicon_Positions'][0], row['Amplicon_Positions'][1], args.segment
        ) if row['Amplicon_Positions'] else 0.0, 
        axis=1
    )
    df['Amplicon_Dropout_Percentage'] = df.apply(
        lambda row: calculate_amplicon_dropout(
            row['Amplicon_Positions'][0], row['Amplicon_Positions'][1], args.segment
        ) if row['Amplicon_Positions'] else 100.0,
        axis=1
    )
    
    # Check primer quality
    print("Evaluating primer quality...")
    df['Quality_Flag'] = df.apply(check_primer_quality, axis=1)
    
    # Plot primers to target sequence
    print("Generating primer map visualizations...")
    plot_output = os.path.join(args.output, f"{args.segment}_primer_map.png")
    plot_primers_to_target(df, consensus_seq, plot_output, args.segment)
    
    # Save results to CSV
    output_csv = os.path.join(args.output, f"{args.segment}_primer_evaluation.csv")
    df.to_csv(output_csv, index=False, encoding='utf-8')
    print(f"Results saved to {output_csv}")
    
    # Print summary
    passing_primers = df[df['Quality_Flag'] == 'Pass']
    print(f"\nSummary: {len(passing_primers)} of {len(df)} primer pairs pass all quality checks")
    
    if not passing_primers.empty:
        print("\nPassing primer pairs:")
        for _, row in passing_primers.iterrows():
            print(f"  - {row['F Name']} / {row['Reverse Name']}")
    
    svg_plot_output = plot_output.replace('.png', '.svg')
    amplicons_only_output = plot_output.replace('.png', '_amplicons_only.png')
    amplicons_only_svg = amplicons_only_output.replace('.png', '.svg')
    
    print(f"\nPrimer map visualizations saved to:")
    print(f"  - {plot_output} (PNG)")
    print(f"  - {svg_plot_output} (SVG)")
    print(f"  - {amplicons_only_output} (PNG, amplicons only)")
    print(f"  - {amplicons_only_svg} (SVG, amplicons only)")

if __name__ == "__main__":
    main() 