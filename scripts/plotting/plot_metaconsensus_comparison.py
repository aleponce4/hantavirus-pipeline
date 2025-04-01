#!/usr/bin/env python3
"""
Simple script to create an alignment visualization between original reference and metaconsensus
using MAFFT and pyMSAviz
"""

import os
import sys
import subprocess
from pathlib import Path
from Bio import SeqIO

# Set non-interactive Agg backend for matplotlib
import matplotlib
matplotlib.use('Agg')

# Import pymsaviz
try:
    from pymsaviz import MsaViz
except ImportError:
    print("Installing pymsaviz...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pymsaviz"])
    from pymsaviz import MsaViz

def create_alignment(segment):
    """Create an alignment between reference and metaconsensus using MAFFT"""
    # Setup paths
    segment_dir = Path(f"data/references/{segment}")
    plots_dir = Path("results/plots")
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # Find original reference
    reference_files = list(segment_dir.glob("*.fasta"))
    original_ref = None
    
    for ref in reference_files:
        if "metaconsensus" not in ref.name:
            original_ref = ref
            break
    
    if not original_ref:
        print(f"No original reference found for {segment}")
        return None
    
    # Get metaconsensus
    metaconsensus = segment_dir / "metaconsensus.fasta"
    
    if not metaconsensus.exists():
        print(f"No metaconsensus found for {segment}")
        return None
    
    # Create combined file for alignment
    combined_fasta = plots_dir / f"{segment}_sequences.fasta"
    
    # Write sequences with clear IDs
    with open(combined_fasta, 'w') as f:
        # Original reference
        ref_record = list(SeqIO.parse(original_ref, "fasta"))[0]
        f.write(f">Original_Reference\n{ref_record.seq}\n")
        
        # Metaconsensus
        meta_record = list(SeqIO.parse(metaconsensus, "fasta"))[0]
        f.write(f">Metaconsensus\n{meta_record.seq}\n")
    
    # Create alignment with MAFFT
    alignment_file = plots_dir / f"{segment}_aligned.fasta"
    
    try:
        # Run MAFFT
        cmd = ["mafft", "--auto", "--quiet", "--preservecase", str(combined_fasta)]
        with open(alignment_file, "w") as out:
            subprocess.run(cmd, stdout=out, check=True)
        
        # Clean up intermediate file
        os.remove(combined_fasta)
        
        return alignment_file
    
    except subprocess.CalledProcessError as e:
        print(f"Error running MAFFT: {e}")
        return None
    except FileNotFoundError:
        print("MAFFT not found. Please install MAFFT.")
        return None

def visualize_alignment(alignment_file, segment):
    """Create a visualization of the alignment using pyMSAviz"""
    if alignment_file is None:
        return
    
    plots_dir = Path("results/plots")
    output_file = plots_dir / f"{segment}_alignment.png"
    
    try:
        # Calculate the appropriate wrap length based on sequence length
        # Get sequence length
        seq_length = 0
        with open(alignment_file, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                seq_length = len(record.seq)
                break
        
        # Set wrap length to create multiple rows (aim for 3-5 rows)
        wrap_length = min(100, max(60, seq_length // 4))
        
        # Create the visualization with basic parameters
        # Use only parameters that are definitely supported
        mv = MsaViz(
            alignment_file,
            wrap_length=wrap_length,  # Set wrap length for multiple rows
            show_count=True,          # Show position counts
            show_grid=True,           # Show grid lines
            show_consensus=True       # Show consensus sequence
        )
        
        # Save the figure with higher resolution
        mv.savefig(str(output_file), dpi=300)
        
        print(f"Alignment visualization saved to {output_file}")
        return output_file
    except Exception as e:
        print(f"Error creating visualization: {e}")
        return None

def main():
    """Main function to create alignment visualizations for segments"""
    segments = ["L_segment", "M_segment", "S_segment"]
    
    for segment in segments:
        print(f"\nProcessing {segment}...")
        
        # Create alignment
        alignment_file = create_alignment(segment)
        
        # Visualize if alignment was successful
        if alignment_file:
            visualize_alignment(alignment_file, segment)
        else:
            print(f"Skipping visualization for {segment}")

if __name__ == "__main__":
    main() 