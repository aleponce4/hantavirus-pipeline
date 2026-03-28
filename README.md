# Hantavirus NGS Analysis Workflow

Workflow for processing Hantavirus Illumina sequencing data including quality control, mapping, variant calling, consensus generation, and basic visualization. Developed to standardize analysis of segmented viral genomes (L, M, and S segments) in research studies.

## Overview

This workflow performs:

1. Read preprocessing and trimming
2. Reference-based alignment
3. Variant calling
4. Consensus sequence generation
5. Coverage analysis
6. Primer evaluation

The implementation focuses on reproducible processing of viral sequencing datasets rather than development of new algorithms.

## Workflow Structure

Major components:

• Read QC and adapter trimming  
• Primer trimming  
• Reference alignment  
• Variant calling  
• Consensus generation  
• Coverage analysis  
• Primer evaluation  

## Installation

### Requirements

• Conda or Miniconda  
• Git (optional)  

### Setup

Clone repository:

```bash
git clone https://github.com/aleponce4/hantavirus-pipeline
cd hantavirus_pipeline
````

Run setup:

```bash
./setup.sh
```

The setup script:

• Creates the conda environment
• Installs required dependencies
• Creates directory structure

## Input Data

Expected directory structure:

### FASTQ files

Place paired-end reads in:

```
data/raw_reads/
```

Naming convention:

```
SAMPLE_R1.fastq.gz
SAMPLE_R2.fastq.gz
```

### Reference sequences

Organize references by genome segment:

```
data/references/

L_segment/
M_segment/
S_segment/
```

Each segment should contain:

• FASTA reference sequence
• GFF3 annotation file

### Primer file (optional)

```
data/primers/Primers.csv
```

Required columns:

• Primer name
• Region
• Sequence

Example:

```
Name,Region,Sequence
SF1,1-32F,TAGTAGTAGACTCCTTGAGAAGCTACT
```

## Running the Workflow

Activate environment:

```bash
conda activate De_Novo_pipeline
```

Run pipeline:

```bash
bash run_pipeline.sh
```

Run multiple samples sequentially:

```bash
for sample in sample1 sample2 sample3
do
    bash run_pipeline.sh
done
```

## HPC Usage

The workflow supports execution on multi-core systems and HPC environments.

Features:

• Automatic CPU detection
• Sample-level parallel processing
• Resource-aware thread allocation
• Job tracking

Example SLURM submission:

```bash
#SBATCH --job-name=hantavirus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=24:00:00

module load conda
conda activate De_Novo_pipeline

cd /path/to/hantavirus_pipeline

bash run_pipeline.sh
```

## Processing Strategy

The workflow uses a two-pass alignment strategy.

### First pass

• Mapping to reference
• Variant calling
• Initial consensus generation

### Sample classification

Samples with low average coverage are classified as negative and excluded from second-pass refinement.

### Second pass (positive samples)

• Alignment to metaconsensus
• Read extraction
• Re-alignment to original reference
• Final variant calling
• Consensus comparison

This improves mapping while maintaining consistent genomic coordinates.

## Output Structure

Results are written to:

```
results/

first_pass/
second_pass/
negative_samples/
trimmed/
plots/
primer_evaluation/
```

Description:

• **first_pass/** – Initial processing results
• **second_pass/** – Refined processing results
• **negative_samples/** – Low coverage samples
• **trimmed/** – Processed reads
• **plots/** – Coverage visualizations
• **primer_evaluation/** – Primer analysis results

## Visualization Outputs

Generated outputs include:

• Coverage plots
• Reference vs consensus alignment comparisons
• Primer location maps

## Primer Evaluation

Primer evaluation module analyzes:

• Primer sequence properties
• Binding efficiency to consensus sequences
• Potential problematic primers

Outputs include summary reports and suggested improvements.

## Troubleshooting

Common checks:

• Verify FASTQ naming conventions
• Confirm conda environment is active
• Check logs directory for errors
• Ensure FASTA and GFF files match
• Verify input directory structure


