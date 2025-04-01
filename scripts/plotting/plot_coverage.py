#!/usr/bin/env python3
"""
Script to generate coverage plots for multiple samples in a CLC Workbench-like style
"""

import os
import sys
import subprocess
import argparse
import csv
from pathlib import Path
import numpy as np
import pysam

# Set non-interactive Agg backend for matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

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

def calculate_average_coverage(bam_file, min_quality=0):
    """Calculate the true average coverage (not capped) for a BAM file"""
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            # Get reference length from the first reference
            ref_name = bam.references[0]
            ref_length = bam.lengths[0]
            
            # Generate coverage array
            coverage = np.zeros(ref_length)
            for pileupcol in bam.pileup(ref_name, min_mapping_quality=min_quality):
                coverage[pileupcol.pos] = pileupcol.n
            
            # Calculate average including zero coverage positions
            avg_coverage = np.mean(coverage)
            return avg_coverage
    except Exception as e:
        print(f"Error calculating average coverage for {bam_file}: {e}")
        return 0.0

def find_samples_and_bams(results_dir, segment):
    """Find all sample BAM files for a specific segment"""
    samples_bams = {}
    
    try:
        results_path = Path(results_dir)
        for sample_dir in results_path.iterdir():
            if not sample_dir.is_dir() or sample_dir.name == "plots" or sample_dir.name == "first_pass" or sample_dir.name == "negative_samples":
                continue
            
            # Check for alignment directory
            align_dir = sample_dir / "alignment" / segment
            if not align_dir.exists() or not align_dir.is_dir():
                continue
            
            # Look for BAM file
            bam_file = align_dir / f"{sample_dir.name}.bam"
            if bam_file.exists():
                samples_bams[sample_dir.name] = str(bam_file)
    
    except Exception as e:
        print(f"Error finding samples: {e}")
    
    return samples_bams

def get_primer_positions(segment):
    """Extract primer binding positions from the Primers.csv file"""
    primer_pairs = []
    primer_file = Path("data/primers/Primers.csv")
    
    if not primer_file.exists():
        print(f"Warning: Primers file not found at {primer_file}")
        return primer_pairs
    
    try:
        with open(primer_file, 'r') as f:
            reader = csv.reader(f)
            next(reader)  # Skip header
            
            prefix = segment[0].upper()  # Get S, M, or L from segment name
            
            for row in reader:
                if len(row) < 10:  # Skip incomplete rows
                    continue
                    
                f_name = row[0]
                f_region = row[1]
                r_name = row[7]
                r_region = row[8]
                
                # Only process primers for the current segment
                if not (f_name.startswith(prefix + "F") or r_name.startswith(prefix + "R")):
                    continue
                
                forward_pos = None
                reverse_pos = None
                
                # Process forward primer position
                if f_region and "-" in f_region:
                    try:
                        # Extract start-end from "1-32F" format
                        pos = f_region.split('-')[0]
                        if pos.isdigit():
                            forward_pos = int(pos)
                    except Exception as e:
                        print(f"Warning: Could not parse forward primer position from {f_region}: {e}")
                
                # Process reverse primer position
                if r_region and "-" in r_region:
                    try:
                        # Extract start-end from "452-480R" format
                        pos = r_region.split('-')[0]
                        if pos.isdigit():
                            reverse_pos = int(pos)
                    except Exception as e:
                        print(f"Warning: Could not parse reverse primer position from {r_region}: {e}")
                
                # Only add if we have both positions
                if forward_pos is not None and reverse_pos is not None:
                    primer_pairs.append({
                        'forward_pos': forward_pos,
                        'reverse_pos': reverse_pos,
                        'forward_name': f_name,
                        'reverse_name': r_name
                    })
    
    except Exception as e:
        print(f"Error reading primers file: {e}")
    
    return primer_pairs

def plot_coverage(samples_bams, output_file, segment, width=10, height_per_sample=2, 
                 min_quality=0, log_scale=False, low_coverage_threshold=10,
                 max_coverage=1000):
    """Generate a multi-sample coverage plot in CLC Workbench style"""
    n_samples = len(samples_bams)
    if n_samples == 0:
        print(f"No samples found for {segment}")
        return False
    
    # Get primer positions for this segment
    primer_pairs = get_primer_positions(segment)
    
    # Define a colormap for primer pairs
    color_cycle = plt.cm.tab20(np.linspace(0, 1, max(10, len(primer_pairs))))
    
    # Set up the figure
    fig, axes = plt.subplots(nrows=n_samples, ncols=1, 
                             figsize=(width, height_per_sample * n_samples),
                             sharex=True)
    
    # Make sure axes is always a list even with one sample
    if n_samples == 1:
        axes = [axes]
    
    # Set up Matplotlib style for clean look
    plt.style.use('seaborn-v0_8-whitegrid')
    
    # Process each sample
    for i, (sample, bam_file) in enumerate(samples_bams.items()):
        # Get coverage data - capped at max_coverage
        coverage = get_coverage_from_bam(bam_file, min_quality, max_coverage)
        if len(coverage) == 0:
            continue
        
        # Calculate true average coverage for display in title (not capped)
        avg_coverage = calculate_average_coverage(bam_file, min_quality)
        
        # Plot coverage
        ax = axes[i]
        x = np.arange(1, len(coverage) + 1)  # 1-based positions
        
        # Fill area under the curve with blue
        ax.fill_between(x, 0, coverage, color='#3366CC', alpha=0.6)
        
        # Plot the coverage line
        ax.plot(x, coverage, '-', color='#1A3C6E', linewidth=1.0)
        
        # Add low coverage highlighting
        low_cov_mask = coverage < low_coverage_threshold
        if np.any(low_cov_mask):
            ax.fill_between(x, 0, coverage, where=low_cov_mask, color='#FF9999', alpha=0.8, 
                           label=f'<{low_coverage_threshold}x')
        
        # Add primer position markers colored by primer pair
        legend_handles = []
        for idx, pair in enumerate(primer_pairs):
            color = color_cycle[idx % len(color_cycle)]
            
            # Draw forward primer
            if pair['forward_pos'] <= len(coverage):
                # Draw forward primer positions as triangles at the top of the plot (using pair color)
                ax.plot(pair['forward_pos'], max_coverage, marker='^', color=color, markersize=8, 
                       markeredgecolor='black', alpha=0.8)
                
                # Add subtle vertical line
                ax.axvline(x=pair['forward_pos'], color=color, linestyle=':', alpha=0.3, linewidth=0.8)
            
            # Draw reverse primer
            if pair['reverse_pos'] <= len(coverage):
                # Draw reverse primer positions as triangles pointing down at the top of the plot (using pair color)
                ax.plot(pair['reverse_pos'], max_coverage, marker='v', color=color, markersize=8, 
                       markeredgecolor='black', alpha=0.8)
                
                # Add subtle vertical line
                ax.axvline(x=pair['reverse_pos'], color=color, linestyle=':', alpha=0.3, linewidth=0.8)
            
            # Add to legend if both primers are within range
            if pair['forward_pos'] <= len(coverage) and pair['reverse_pos'] <= len(coverage):
                # Add a label with F and R primer names
                label = f"{pair['forward_name']}-{pair['reverse_name']}"
                handle = plt.Line2D([0], [0], marker='o', color=color, markeredgecolor='black',
                                   markersize=6, label=label, linestyle='-', linewidth=1)
                legend_handles.append(handle)
        
        # Set title and labels with average coverage
        ax.set_title(f"{sample} [Avg: {avg_coverage:.1f}x]", fontsize=10, loc='left', fontweight='bold')
        
        # Only show y-label for the middle plot
        if i == n_samples // 2:
            ax.set_ylabel('Coverage', fontsize=10)
        
        # Set y-axis to log scale if requested
        if log_scale:
            ax.set_yscale('log')
        
        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Add light grid lines
        ax.grid(True, linestyle='--', alpha=0.7, which='major')
        
        # Add a subtle horizontal line at the low coverage threshold
        ax.axhline(y=low_coverage_threshold, linestyle='--', color='#FF0000', alpha=0.5, 
                  linewidth=0.8)
        
        # Set y-ticks and limits
        if log_scale:
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
            ax.set_ylim([0.5, max_coverage * 1.1])
        else:
            # For linear scale, set max y limit to max_coverage
            ax.set_ylim([0, max_coverage * 1.1])
            ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
    # Set the x-axis label on the bottom-most subplot
    axes[-1].set_xlabel('Position', fontsize=10)
    
    # Add legend for primer pairs if we have any
    if legend_handles:
        # If we have more than 10 pairs, it's better to place the legend outside
        if len(legend_handles) > 10:
            # Place the legend to the right of the figure
            fig.subplots_adjust(right=0.8)
            legend = fig.legend(handles=legend_handles, loc='center right', 
                               bbox_to_anchor=(0.98, 0.5), title="Primer Pairs",
                               frameon=True, framealpha=0.7, fontsize=8)
        else:
            # For fewer pairs, place it in the top-right of the first subplot
            legend = axes[0].legend(handles=legend_handles, loc='upper right', 
                                  title="Primer Pairs",
                                  frameon=True, framealpha=0.7, fontsize=8)
    
    # Set the figure title
    plt.suptitle(f"{segment}", fontsize=14, fontweight='bold')
    
    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(top=0.95, hspace=0.3)
    
    # Save the figure
    try:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Coverage plot saved to: {output_file}")
        plt.close()
        return True
    except Exception as e:
        print(f"Error saving plot: {e}")
        return False

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Generate coverage plots for samples")
    parser.add_argument("--segment", default="S_segment", help="Segment to plot (default: S_segment)")
    parser.add_argument("--results_dir", default="results", help="Directory containing results (default: results)")
    parser.add_argument("--width", type=float, default=10, help="Width of the plot in inches (default: 10)")
    parser.add_argument("--height_per_sample", type=float, default=2, 
                      help="Height per sample in inches (default: 2)")
    parser.add_argument("--min_quality", type=int, default=0, 
                      help="Minimum mapping quality to include (default: 0)")
    parser.add_argument("--low_coverage", type=int, default=10, 
                      help="Threshold for low coverage highlighting (default: 10)")
    parser.add_argument("--max_coverage", type=int, default=1000,
                      help="Max coverage cap (default: 1000)")
    parser.add_argument("--log_scale", action="store_true", help="Use log scale for y-axis")
    parser.add_argument("--output", help="Output file path (default: auto-generated)")
    
    # Check if the script was called without args (meaning from pipeline)
    if len(sys.argv) == 1:
        # Find available segments
        results_dir = "results"
        segments = []
        
        for segment in ["L_segment", "M_segment", "S_segment"]:
            samples_bams = find_samples_and_bams(results_dir, segment)
            if samples_bams:
                segments.append(segment)
                
                # Create output directory
                plots_dir = Path("results/plots")
                plots_dir.mkdir(parents=True, exist_ok=True)
                
                # Generate plot
                output_file = str(plots_dir / f"{segment}_coverage.png")
                plot_coverage(samples_bams, output_file, segment)
        
        if not segments:
            print("No segments with BAM files found")
    else:
        # Parse arguments
        args = parser.parse_args()
        
        # Find samples and their BAM files
        samples_bams = find_samples_and_bams(args.results_dir, args.segment)
        
        if not samples_bams:
            print(f"No samples found with BAM files for {args.segment}")
            sys.exit(1)
        
        # Set output file if not specified
        if not args.output:
            plots_dir = Path("results/plots")
            plots_dir.mkdir(parents=True, exist_ok=True)
            output_file = str(plots_dir / f"{args.segment}_coverage.png")
        else:
            output_file = args.output
        
        # Generate plot
        success = plot_coverage(
            samples_bams, 
            output_file, 
            args.segment,
            width=args.width,
            height_per_sample=args.height_per_sample,
            min_quality=args.min_quality,
            log_scale=args.log_scale,
            low_coverage_threshold=args.low_coverage,
            max_coverage=args.max_coverage
        )
        
        if not success:
            sys.exit(1)

if __name__ == "__main__":
    main() 