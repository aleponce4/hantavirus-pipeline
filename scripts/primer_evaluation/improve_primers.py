#!/usr/bin/env python3
"""
Primer Improvement Script for Hantavirus Pipeline

This script analyzes flagged primers from the evaluation report and suggests improvements:
1. First fixes binding site mismatches using the reference sequence
2. Then attempts to address other issues like Tm, GC content, etc.
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt
import seqfold
import re

# Import functions from evaluate_primers.py
try:
    from evaluate_primers import calculate_primer_properties, predict_gibbs_free_energy, verify_primer_binding_site, run_blast_and_get_evalue
except ImportError:
    # Handle case when script is run from a different directory
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from evaluate_primers import calculate_primer_properties, predict_gibbs_free_energy, verify_primer_binding_site, run_blast_and_get_evalue

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Suggest improvements for flagged primers')
    parser.add_argument('--eval_report', required=True, help='Path to primer evaluation CSV report')
    parser.add_argument('--reference', required=True, help='Path to reference sequence FASTA')
    parser.add_argument('--output', required=True, help='Path to output CSV with suggested improvements')
    parser.add_argument('--blast_db', help='Path to BLAST database for primer specificity checks')
    return parser.parse_args()

def get_perfect_binding_primer(reference_seq, position, is_reverse=False):
    """
    Extract the exact sequence from reference at the intended binding site.
    
    Args:
        reference_seq (str): Reference sequence
        position (tuple): (start, end) in 0-based coordinates
        is_reverse (bool): Whether this is a reverse primer
        
    Returns:
        str: Perfect-binding primer sequence
    """
    start, end = position
    
    # Ensure positions are within bounds
    start = max(0, start)
    end = min(len(reference_seq) - 1, end)
    
    # Extract the sequence
    primer_seq = reference_seq[start:end+1]
    
    # For reverse primers, we need to return the reverse complement
    if is_reverse:
        primer_seq = str(Seq(primer_seq).reverse_complement())
        
    return primer_seq

def check_primer_issues(primer_seq, binding_match, tm, gc_content, gibbs_energy, blast_evalue=None):
    """
    Check for issues with a primer and return a list of problems.
    
    Args:
        primer_seq (str): Primer sequence
        binding_match (float): Binding match percentage
        tm (float): Melting temperature
        gc_content (float): GC content
        gibbs_energy (float): Gibbs free energy
        blast_evalue (float, optional): BLAST E-value for specificity check
        
    Returns:
        list: List of issues with the primer
    """
    issues = []
    
    # Check binding
    if binding_match < 100:
        issues.append(f"Binding match ({binding_match:.1f}%)")
    
    # Check GC content
    if not (0.35 <= gc_content <= 0.65):
        if gc_content < 0.35:
            issues.append(f"Low GC content ({gc_content:.2f})")
        else:
            issues.append(f"High GC content ({gc_content:.2f})")
    
    # Check melting temperature - using tighter range of 60-68°C
    if not (60 <= tm <= 68):
        if tm < 60:
            issues.append(f"Low Tm ({tm:.1f}°C)")
        else:
            issues.append(f"High Tm ({tm:.1f}°C)")
    
    # Check Gibbs Free Energy
    if gibbs_energy < -6:
        issues.append(f"Hairpin potential (ΔG={gibbs_energy:.2f})")
    
    # Check primer length - Updated to use 21 as the minimum
    if not (21 <= len(primer_seq) <= 35):
        if len(primer_seq) < 21:
            issues.append(f"Too short ({len(primer_seq)}bp)")
        else:
            issues.append(f"Too long ({len(primer_seq)}bp)")
    
    # Check for runs of the same base
    for base in 'ACGT':
        if base * 5 in primer_seq:
            issues.append(f"Run of {base}s")
            break
    
    # Check BLAST specificity if provided
    if blast_evalue is not None and blast_evalue > 1e-4:
        issues.append(f"Poor specificity (E-value={blast_evalue:.1e})")
    
    return issues

def improve_primer_for_binding(original_primer, reference_seq, position, is_reverse):
    """
    Replace a primer with perfect binding sequence from reference.
    
    Args:
        original_primer (str): Original primer sequence
        reference_seq (str): Reference sequence
        position (tuple): (start, end) in 0-based coordinates
        is_reverse (bool): Whether this is a reverse primer
        
    Returns:
        dict: Improved primer information
    """
    # Get perfect binding sequence
    perfect_primer = get_perfect_binding_primer(reference_seq, position, is_reverse)
    
    # Clean up any whitespace or non-standard characters
    original_primer = ''.join(c for c in original_primer if c.upper() in 'ACGT')
    
    # Calculate number of changes
    changes = sum(a != b for a, b in zip(original_primer, perfect_primer))
    
    # Calculate new properties
    gc_content, tm = calculate_primer_properties(perfect_primer)
    gibbs_energy = predict_gibbs_free_energy(perfect_primer)
    
    return {
        'sequence': perfect_primer,
        'changes': changes,
        'gc_content': gc_content,
        'tm': tm,
        'gibbs_energy': gibbs_energy,
        'notes': f"Modified for perfect binding ({changes} changes)"
    }

def optimize_primer_tm(primer_seq, reference_seq, position, is_reverse, target_tm=65):
    """
    Optimize primer melting temperature using a combination of:
    1. Trimming to ideal length (25bp)
    2. Sliding the binding position slightly (+/- 3bp)
    
    This ensures primers will have optimal Tm in the 60-68°C range.
    
    Args:
        primer_seq (str): Current primer sequence
        reference_seq (str): Reference sequence
        position (tuple): Current binding position (start, end)
        is_reverse (bool): Whether this is a reverse primer
        target_tm (float): Target Tm (default 65°C)
        
    Returns:
        tuple: (optimized_primer, new_tm, new_gc, changes_made)
    """
    # Clean the primer sequence
    primer_seq = ''.join(c for c in primer_seq if c.upper() in 'ACGT')
    
    # Check original properties
    orig_gc, orig_tm = calculate_primer_properties(primer_seq)
    
    # If already in optimal range, return original
    if 60 <= orig_tm <= 68:
        return primer_seq, orig_tm, orig_gc, 0
    
    # Strategy 1: Trim to ideal 25bp length if longer
    if len(primer_seq) > 25:
        if is_reverse:
            # For reverse primers, keep the 3' end (right side)
            trimmed = primer_seq[-25:]
        else:
            # For forward primers, keep the 3' end (right side) 
            trimmed = primer_seq[-25:] if len(primer_seq) <= 30 else primer_seq[:25]
        
        gc, tm = calculate_primer_properties(trimmed)
        
        # If trimming fixed the issue, return it
        if 60 <= tm <= 68:
            return trimmed, tm, gc, len(primer_seq) - len(trimmed)
    else:
        # Use original if not longer than ideal
        trimmed = primer_seq
        gc, tm = orig_gc, orig_tm
    
    # Strategy 2: If we have a reference position, try sliding the binding site
    if position and len(reference_seq) > 0:
        start, end = position
        
        # Prepare to find the best candidate
        candidates = []
        
        # Try the standard window with offsets: -3, +3
        for offset in [-3, 0, 3]:
            new_start = max(0, start + offset)
            new_end = min(len(reference_seq) - 1, end + offset)
            
            # Skip if invalid range
            if new_end - new_start + 1 < 21:
                continue
                
            # Extract sequence from reference at this position
            if is_reverse:
                ref_seq = str(Seq(reference_seq[new_start:new_end+1]).reverse_complement())
            else:
                ref_seq = reference_seq[new_start:new_end+1]
            
            # If longer than 25bp, trim to ideal length
            if len(ref_seq) > 25:
                ref_seq = ref_seq[-25:] if is_reverse else ref_seq[:25]
            
            # Calculate properties and store as candidate
            cand_gc, cand_tm = calculate_primer_properties(ref_seq)
            tm_diff = abs(cand_tm - target_tm)
            
            candidates.append((ref_seq, cand_tm, cand_gc, tm_diff))
        
        # Find the best candidate based on Tm closest to target
        if candidates:
            candidates.sort(key=lambda x: x[3])  # Sort by Tm difference
            best_primer, best_tm, best_gc, _ = candidates[0]
            
            # If best candidate has better Tm than trimmed, use it
            if abs(best_tm - target_tm) < abs(tm - target_tm):
                changes = sum(1 for a, b in zip(primer_seq, best_primer) if a != b)
                if len(primer_seq) != len(best_primer):
                    changes += abs(len(primer_seq) - len(best_primer))
                return best_primer, best_tm, best_gc, changes
    
    # Default to the trimmed version if sliding didn't work
    changes = len(primer_seq) - len(trimmed)
    return trimmed, tm, gc, changes

def improve_primer(original_primer, binding_match, tm, gc_content, gibbs_energy, 
                   reference_seq, position, is_reverse, blast_db=None, max_iterations=3):
    """
    Attempt to improve a primer based on its issues.
    
    Args:
        original_primer (str): Original primer sequence
        binding_match (float): Binding match percentage
        tm (float): Melting temperature
        gc_content (float): GC content
        gibbs_energy (float): Gibbs free energy
        reference_seq (str): Reference sequence
        position (tuple): (start, end) in 0-based coordinates
        is_reverse (bool): Whether this is a reverse primer
        blast_db (str, optional): Path to BLAST database for specificity check
        max_iterations (int): Maximum number of improvement iterations
        
    Returns:
        dict: Improved primer information
    """
    # Clean up original primer - remove spaces and non-ACGT characters
    original_primer = ''.join(c for c in original_primer if c.upper() in 'ACGT')
    
    # Get original BLAST E-value if database provided
    original_evalue = None
    if blast_db:
        original_evalue = run_blast_and_get_evalue(original_primer, blast_db)
    
    # Check for initial issues
    initial_issues = check_primer_issues(original_primer, binding_match, tm, gc_content, gibbs_energy, original_evalue)
    
    # If no issues, return the original
    if not initial_issues:
        return {
            'sequence': original_primer,
            'changes': 0,
            'notes': "No changes needed - primer passes all checks"
        }
    
    # If binding match is less than 100%, replace with perfect binding sequence
    if binding_match < 100 and position is not None:
        improved = improve_primer_for_binding(original_primer, reference_seq, position, is_reverse)
        
        # Check BLAST E-value for improved primer
        improved_evalue = None
        if blast_db:
            improved_evalue = run_blast_and_get_evalue(improved['sequence'], blast_db)
            improved['blast_evalue'] = improved_evalue
        
        # Check if the perfect binding primer has no other issues
        new_issues = check_primer_issues(
            improved['sequence'], 
            100, # Perfect binding 
            improved['tm'], 
            improved['gc_content'], 
            improved['gibbs_energy'],
            improved_evalue
        )
        
        # If no other issues, return the perfect binding primer
        if not new_issues:
            improved['notes'] = "Fixed binding issue - primer now passes all checks"
            return improved
        
        # Otherwise, continue improving from this base
        current_primer = improved['sequence']
        current_gc = improved['gc_content']
        current_tm = improved['tm']
        current_energy = improved['gibbs_energy']
        current_evalue = improved_evalue
        changes = improved['changes']
        notes = "Fixed binding issue, but still has other problems: " + ", ".join(new_issues)
    else:
        # Start with original if binding is already good
        current_primer = original_primer
        current_gc = gc_content
        current_tm = tm
        current_energy = gibbs_energy
        current_evalue = original_evalue
        changes = 0
        notes = "Binding is good, addressing other issues: " + ", ".join(initial_issues)
    
    # Store best version so far
    best_primer = {
        'sequence': current_primer,
        'gc_content': current_gc,
        'tm': current_tm,
        'gibbs_energy': current_energy,
        'blast_evalue': current_evalue,
        'issues': check_primer_issues(current_primer, 100, current_tm, current_gc, current_energy, current_evalue),
        'changes': changes
    }
    
    # Improvement iterations
    for iteration in range(max_iterations):
        # Get current issues
        current_issues = check_primer_issues(current_primer, 100, current_tm, current_gc, current_energy, current_evalue)
        
        # If no issues, we're done
        if not current_issues:
            notes += f" - Fixed after {iteration+1} iterations"
            return {
                'sequence': current_primer,
                'changes': changes,
                'gc_content': current_gc,
                'tm': current_tm,
                'gibbs_energy': current_energy,
                'blast_evalue': current_evalue,
                'notes': notes
            }
        
        # Try to fix each issue
        improved = False
        
        # Check for temperature issues (high or low Tm)
        if any("Tm" in issue for issue in current_issues):
            # Use the new optimize_primer_tm function to address Tm issues
            optimized_primer, new_tm, new_gc, tm_changes = optimize_primer_tm(
                current_primer, reference_seq, position, is_reverse
            )
            
            if optimized_primer != current_primer:
                current_primer = optimized_primer
                current_gc = new_gc
                current_tm = new_tm
                current_energy = predict_gibbs_free_energy(current_primer)
                if blast_db:
                    current_evalue = run_blast_and_get_evalue(current_primer, blast_db)
                changes += tm_changes
                improved = True
                print(f"  Optimized primer for Tm: {current_tm:.2f}°C, length now {len(current_primer)}bp")
        
        # Handle other issues if Tm wasn't fixed or wasn't an issue
        if not improved:
            # Handle specificity issues if BLAST DB provided
            if blast_db and any("Poor specificity" in issue for issue in current_issues):
                # Try extending the primer to improve specificity (if not too long already)
                if len(current_primer) < 30 and position is not None:
                    # Get the full region from reference
                    full_seq = get_perfect_binding_primer(reference_seq, position, is_reverse)
                    if len(full_seq) > len(current_primer):
                        # For forward primers, try adding bases to 3' end
                        # For reverse primers, try adding bases to 3' end (biologically relevant end)
                        if is_reverse:
                            # Add to 3' end (right) of the primer from the full sequence
                            extended = current_primer + full_seq[len(current_primer):min(len(current_primer)+3, len(full_seq))]
                        else:
                            # Add to 3' end (right) of the primer from the full sequence
                            extended = current_primer + full_seq[len(current_primer):min(len(current_primer)+3, len(full_seq))]
                        
                        # Check if extended primer has better specificity
                        extended_evalue = run_blast_and_get_evalue(extended, blast_db)
                        if extended_evalue < current_evalue:
                            current_primer = extended
                            current_gc, current_tm = calculate_primer_properties(current_primer)
                            current_energy = predict_gibbs_free_energy(current_primer)
                            current_evalue = extended_evalue
                            changes += 1
                            improved = True
                            print(f"  Extended primer to improve specificity: E-value now {current_evalue:.1e}")
            
            # Then GC content
            elif any("GC content" in issue for issue in current_issues):
                target_gc = 0.5  # Midpoint of acceptable range
                new_primer = improve_gc_content(current_primer, target_gc, current_gc)
                if new_primer != current_primer:
                    current_primer = new_primer
                    current_gc, current_tm = calculate_primer_properties(current_primer)
                    current_energy = predict_gibbs_free_energy(current_primer)
                    if blast_db:
                        current_evalue = run_blast_and_get_evalue(current_primer, blast_db)
                    changes += 1
                    improved = True
            
            # Then hairpins
            elif any("Hairpin" in issue for issue in current_issues):
                new_primer = fix_hairpins(current_primer, current_energy)
                if new_primer != current_primer:
                    current_primer = new_primer
                    current_gc, current_tm = calculate_primer_properties(current_primer)
                    current_energy = predict_gibbs_free_energy(current_primer)
                    if blast_db:
                        current_evalue = run_blast_and_get_evalue(current_primer, blast_db)
                    changes += 1
                    improved = True
            
            # Then length issues
            elif any(("Too short" in issue or "Too long" in issue) for issue in current_issues):
                target_length = 25  # Ideal length
                new_primer = adjust_length(current_primer, target_length, is_reverse)
                if new_primer != current_primer:
                    current_primer = new_primer
                    current_gc, current_tm = calculate_primer_properties(current_primer)
                    current_energy = predict_gibbs_free_energy(current_primer)
                    if blast_db:
                        current_evalue = run_blast_and_get_evalue(current_primer, blast_db)
                    changes += 1
                    improved = True
        
        # If we couldn't improve anything, break the loop
        if not improved:
            break
        
        # Update best primer if this one has fewer issues
        current_issue_count = len(check_primer_issues(
            current_primer, 100, current_tm, current_gc, current_energy, current_evalue
        ))
        best_issue_count = len(best_primer['issues'])
        
        if current_issue_count < best_issue_count:
            best_primer = {
                'sequence': current_primer,
                'gc_content': current_gc,
                'tm': current_tm,
                'gibbs_energy': current_energy,
                'blast_evalue': current_evalue,
                'issues': check_primer_issues(current_primer, 100, current_tm, current_gc, current_energy, current_evalue),
                'changes': changes
            }
    
    # Return best primer found
    if best_primer['issues']:
        remaining_issues = ", ".join(best_primer['issues'])
        notes += f" - Made {best_primer['changes']} changes but still has issues: {remaining_issues}"
    else:
        notes += f" - Fixed after {best_primer['changes']} changes"
    
    return {
        'sequence': best_primer['sequence'],
        'changes': best_primer['changes'],
        'gc_content': best_primer.get('gc_content', gc_content),
        'tm': best_primer.get('tm', tm),
        'gibbs_energy': best_primer.get('gibbs_energy', gibbs_energy),
        'blast_evalue': best_primer.get('blast_evalue', None),
        'notes': notes
    }

def improve_gc_content(primer, target_gc, current_gc):
    """
    Modify primer to improve its GC content.
    
    Args:
        primer (str): Current primer sequence
        target_gc (float): Target GC content (midpoint of acceptable range)
        current_gc (float): Current GC content
        
    Returns:
        str: Improved primer sequence
    """
    modified_primer = primer
    
    # Determine if we need to increase or decrease GC
    if current_gc < target_gc:
        # To increase GC, replace A/T with G/C at non-critical positions
        # Avoid first and last 5 bases if possible
        safe_range = range(5, len(primer)-5) if len(primer) > 15 else range(1, len(primer)-1)
        for i in safe_range:
            if primer[i] in 'AT':
                candidate = primer[:i] + ('G' if primer[i] == 'A' else 'C') + primer[i+1:]
                new_gc = calculate_primer_properties(candidate)[0]
                if abs(new_gc - target_gc) < abs(current_gc - target_gc):
                    modified_primer = candidate
                    break
    else:
        # To decrease GC, replace G/C with A/T
        safe_range = range(5, len(primer)-5) if len(primer) > 15 else range(1, len(primer)-1)
        for i in safe_range:
            if primer[i] in 'GC':
                candidate = primer[:i] + ('A' if primer[i] == 'G' else 'T') + primer[i+1:]
                new_gc = calculate_primer_properties(candidate)[0]
                if abs(new_gc - target_gc) < abs(current_gc - target_gc):
                    modified_primer = candidate
                    break
    
    return modified_primer

def fix_hairpins(primer, current_energy):
    """
    Attempt to fix hairpin issues in the primer.
    
    Args:
        primer (str): Current primer sequence
        current_energy (float): Current Gibbs free energy
        
    Returns:
        str: Improved primer sequence
    """
    # Find potential hairpin regions using seqfold
    try:
        # Use a different approach to analyze hairpins that doesn't rely on struct_type
        # This avoids the "Struct object is not subscriptable" error
        fold_data = seqfold.fold(primer)
        
        # If no folding data, return original
        if not fold_data:
            return primer
        
        # Try different positions in the sequence to fix hairpins
        for i in range(len(primer)):
            # Skip the first and last few bases which are critical for binding
            if i < 3 or i > len(primer) - 4:
                continue
                
            base = primer[i]
            alternatives = [b for b in 'ACGT' if b != base]
            
            # Try each alternative base
            for alt_base in alternatives:
                candidate = primer[:i] + alt_base + primer[i+1:]
                new_energy = predict_gibbs_free_energy(candidate)
                
                # If we improved the energy, return this candidate
                if new_energy > current_energy:
                    return candidate
        
    except Exception as e:
        print(f"Error analyzing hairpins: {e}")
    
    # If we couldn't fix it or there was an error, return the original
    return primer

def adjust_length(primer, target_length, is_reverse):
    """
    Adjust primer length to be within acceptable range.
    
    Args:
        primer (str): Current primer sequence
        target_length (int): Target length (ideal length)
        is_reverse (bool): Whether this is a reverse primer
        
    Returns:
        str: Improved primer sequence
    """
    if len(primer) < 21:  # Updated to minimum 21bp
        # Too short - extend the non-critical end
        if is_reverse:
            # For reverse primers, extend 5' end
            # Use a balanced extension
            extension = 'GCAT' * ((21 - len(primer) + 3) // 4)
            return extension[:21-len(primer)] + primer
        else:
            # For forward primers, extend 5' end
            extension = 'GCAT' * ((21 - len(primer) + 3) // 4)
            return extension[:21-len(primer)] + primer
    
    elif len(primer) > 35:
        # Too long - trim the non-critical end
        if is_reverse:
            # For reverse primers, trim 5' end
            return primer[len(primer)-35:]
        else:
            # For forward primers, trim 5' end
            return primer[len(primer)-35:]
    
    # Already within range
    return primer

def extract_position_from_region(region, reference_length):
    """
    Extract start and end positions from region string.
    
    Args:
        region (str): Region string (e.g., "1-32F")
        reference_length (int): Length of reference sequence
        
    Returns:
        tuple: (start, end) in 0-based coordinates, or None if parsing fails
    """
    try:
        if not region or '-' not in region:
            return None
        
        parts = region.split('-')
        start = int(parts[0])
        # Extract digits from second part (removing F/R suffix)
        end = int(''.join(c for c in parts[1] if c.isdigit()))
        
        # Convert to 0-based
        start = start - 1
        end = end - 1
        
        # Validate bounds
        if start < 0 or end >= reference_length or start > end:
            return None
        
        return (start, end)
    except (ValueError, IndexError):
        return None

def suggest_primer_improvements(eval_report_path, reference_path, output_path, blast_db_path=None):
    """
    Process evaluation report and suggest primer improvements.
    
    Args:
        eval_report_path (str): Path to primer evaluation CSV
        reference_path (str): Path to reference sequence FASTA
        output_path (str): Path to output CSV for improved primers
        blast_db_path (str, optional): Path to BLAST database for specificity checks
    """
    # Load evaluation report
    try:
        df = pd.read_csv(eval_report_path)
    except Exception as e:
        print(f"Error reading evaluation report: {e}")
        return False
    
    # Filter to flagged primers
    flagged_df = df[df['Quality_Flag'] != 'Pass'].copy()
    
    if flagged_df.empty:
        print("No flagged primers found in the evaluation report.")
        # Create empty output file with headers
        pd.DataFrame(columns=[
            'Primer_Name', 'Primer_Type', 'Original_Sequence', 'Improved_Sequence',
            'Original_GC', 'Improved_GC', 'Original_Tm', 'Improved_Tm',
            'Original_Gibbs', 'Improved_Gibbs', 'Original_Binding', 'Improved_Binding',
            'Original_E_Value', 'Improved_E_Value', 'Changes', 'Length_Change'
        ]).to_csv(output_path, index=False)
        return True
    
    # Load reference sequence
    try:
        reference_seq = str(SeqIO.read(reference_path, "fasta").seq)
    except Exception as e:
        print(f"Error reading reference sequence: {e}")
        return False
    
    # Use provided BLAST database if available
    blast_db = blast_db_path
    
    # If not provided directly, try to find it automatically
    if not blast_db:
        blast_db_base = reference_path.rsplit('/', 1)[0]
        possible_paths = [
            f"{blast_db_base}/blast_db/{reference_path.rsplit('/', 1)[1].split('.')[0]}_blast_db",
            f"data/blast_db/{reference_path.split('/')[-1].split('.')[0]}_blast_db"
        ]
        
        for path in possible_paths:
            if os.path.exists(path + ".nin") or os.path.exists(path + ".nhr"):
                blast_db = path
                print(f"Found BLAST database: {blast_db}")
                break
    
    if not blast_db:
        print("Warning: No BLAST database found, skipping specificity checks")
    else:
        print(f"Using BLAST database: {blast_db}")
        # Verify BLAST database is accessible
        try:
            test_seq = "ACGTACGTACGT"
            test_result = run_blast_and_get_evalue(test_seq, blast_db)
            print(f"BLAST database test successful: E-value for test sequence = {test_result}")
        except Exception as e:
            print(f"Warning: Error accessing BLAST database: {e}")
            blast_db = None
    
    # Prepare results dataframe
    results = []
    
    # Process each flagged primer pair
    for idx, row in flagged_df.iterrows():
        # Get forward primer info
        f_name = row.get('F Name', '')
        f_seq = row.get('F 5´-3´ Sequence', '')
        if 'F 5´-3´ Sequence' not in row and 'F 5\'-3\' Sequence' in row:
            f_seq = row.get('F 5\'-3\' Sequence', '')
        f_region = row.get('F Region', '')
        f_binding = row.get('F_Binding_Match_Percent', 0)
        f_tm = row.get('F_Melting_Temperature', 0)
        f_gc = row.get('F_GC_Content', 0)
        f_energy = row.get('F_Gibbs_Free_Energy', 0)
        f_evalue = row.get('F_BLAST_E_Value', None)
        
        # Get reverse primer info
        r_name = row.get('Reverse Name', '')
        r_seq = row.get('R 5´-3´ Sequence', '')
        if 'R 5´-3´ Sequence' not in row and 'R 5\'-3\' Sequence' in row:
            r_seq = row.get('R 5\'-3\' Sequence', '')
        r_region = row.get('R Region', '')
        r_binding = row.get('R_Binding_Match_Percent', 0)
        r_tm = row.get('R_Melting_Temperature', 0)
        r_gc = row.get('R_GC_Content', 0)
        r_energy = row.get('R_Gibbs_Free_Energy', 0)
        r_evalue = row.get('R_BLAST_E_Value', None)
        
        # Extract positions
        f_position = extract_position_from_region(f_region, len(reference_seq))
        r_position = extract_position_from_region(r_region, len(reference_seq))
        
        # Improve forward primer
        f_improved = improve_primer(
            f_seq, f_binding, f_tm, f_gc, f_energy,
            reference_seq, f_position, is_reverse=False, blast_db=blast_db
        )
        
        # Improve reverse primer
        r_improved = improve_primer(
            r_seq, r_binding, r_tm, r_gc, r_energy,
            reference_seq, r_position, is_reverse=True, blast_db=blast_db
        )
        
        # Log the improvement notes
        print(f"Forward primer {f_name}: {f_improved['notes']}")
        print(f"Reverse primer {r_name}: {r_improved['notes']}")
        
        # Add forward primer to results (separate row)
        results.append({
            'Primer_Name': f_name,
            'Primer_Type': 'Forward',
            'Original_Sequence': f_seq,
            'Improved_Sequence': f_improved['sequence'],
            'Original_GC': f_gc,
            'Improved_GC': f_improved.get('gc_content', f_gc),
            'Original_Tm': f_tm,
            'Improved_Tm': f_improved.get('tm', f_tm),
            'Original_Gibbs': f_energy,
            'Improved_Gibbs': f_improved.get('gibbs_energy', f_energy),
            'Original_Binding': f_binding,
            'Improved_Binding': 100.0,  # Assuming perfect binding after improvement
            'Original_E_Value': f_evalue,
            'Improved_E_Value': f_improved.get('blast_evalue', None),
            'Changes': f_improved['changes'],
            'Length_Change': len(f_improved['sequence']) - len(f_seq)
        })
        
        # Add reverse primer to results (separate row)
        results.append({
            'Primer_Name': r_name,
            'Primer_Type': 'Reverse',
            'Original_Sequence': r_seq,
            'Improved_Sequence': r_improved['sequence'],
            'Original_GC': r_gc,
            'Improved_GC': r_improved.get('gc_content', r_gc),
            'Original_Tm': r_tm,
            'Improved_Tm': r_improved.get('tm', r_tm),
            'Original_Gibbs': r_energy,
            'Improved_Gibbs': r_improved.get('gibbs_energy', r_energy),
            'Original_Binding': r_binding,
            'Improved_Binding': 100.0,  # Assuming perfect binding after improvement
            'Original_E_Value': r_evalue,
            'Improved_E_Value': r_improved.get('blast_evalue', None),
            'Changes': r_improved['changes'],
            'Length_Change': len(r_improved['sequence']) - len(r_seq)
        })
    
    # Save results
    try:
        pd.DataFrame(results).to_csv(output_path, index=False)
        print(f"Primer improvement suggestions saved to {output_path}")
        print(f"Processed {len(flagged_df)} primer pairs, generated {len(results)} individual primer improvements")
        return True
    except Exception as e:
        print(f"Error saving results: {e}")
        return False

def main():
    args = parse_args()
    
    print(f"Starting primer improvement process...")
    print(f"Evaluation report: {args.eval_report}")
    print(f"Reference sequence: {args.reference}")
    print(f"Output path: {args.output}")
    if args.blast_db:
        print(f"BLAST database: {args.blast_db}")
    
    success = suggest_primer_improvements(
        args.eval_report,
        args.reference,
        args.output,
        args.blast_db
    )
    
    if success:
        print("Primer improvement process completed successfully.")
    else:
        print("Primer improvement process failed.")
        sys.exit(1)

if __name__ == "__main__":
    main() 