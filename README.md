# Hantavirus Pipeline

A pipeline for processing and analyzing Hantavirus NGS data, including quality control, mapping, primer removal, variant calling, and visualization.

## Overview

This pipeline is designed for processing Illumina paired-end sequencing data of Hantavirus genomes. It handles:

- Quality trimming and adapter removal
- Primer trimming
- Reference-based mapping
- Metaconsensus generation
- Variant calling
- Coverage visualization
- Alignment visualization

## Installation

### Prerequisites

- Conda or Miniconda (download from [https://docs.conda.io/projects/miniconda/en/latest/](https://docs.conda.io/projects/miniconda/en/latest/))
- Git (optional, for cloning the repository)

### Setup

1. Clone or download the repository:
   ```bash
   git clone [repository URL]
   cd hantavirus_pipeline
   ```

2. Run the setup script to create the conda environment and necessary directories:
   ```bash
   ./setup.sh
   ```

3. The setup script will:
   - Create a conda environment called `De_Novo_pipeline` with all required dependencies
   - Verify key software installations
   - Create necessary data and results directories

## Usage

### Preparing Input Data

1. Place your paired-end FASTQ files in `data/raw_reads/` using the naming convention:
   - `SAMPLE_R1.fastq.gz` for forward reads
   - `SAMPLE_R2.fastq.gz` for reverse reads

2. Place your reference sequences in `data/references/`:
   - Create a directory for each segment: `L_segment`, `M_segment`, and `S_segment`
   - In each directory, place a FASTA file containing the reference sequence

3. Place your primer information in `data/primers/Primers.csv` (if using primer trimming)

### Running the Pipeline

1. Activate the conda environment:
   ```bash
   conda activate De_Novo_pipeline
   ```

2. Run the pipeline on a sample:
   ```bash
   bash run_pipeline.sh SAMPLE [threads]
   ```
   Where `SAMPLE` is the name of your sample (without the _R1/_R2 suffix) and `threads` is the number of CPU threads to use (default: 4).

3. Alternatively, to run multiple samples sequentially:
   ```bash
   for sample in sample1 sample2 sample3; do
       bash run_pipeline.sh $sample 8
   done
   ```

### Pipeline Steps

The pipeline will execute the following steps:

1. **Preprocessing**: Quality trimming, adapter removal, and primer trimming
2. **Mapping**: Align reads to reference(s) and generate BAM files
3. **Variant calling**: Identify SNPs and indels 
4. **Consensus generation**: Create consensus sequences for each segment
5. **Visualization**: Generate coverage plots and alignment visualizations

### Output Files

Results are organized in the `results/` directory by sample name:

- `results/SAMPLE/trimmed/`: Trimmed and processed FASTQ files
- `results/SAMPLE/alignment/SEGMENT/`: BAM files and statistics
- `results/SAMPLE/variants/SEGMENT/`: Variant calls (VCF)
- `results/SAMPLE/consensus/SEGMENT/`: Consensus sequences
- `results/plots/`: Visualization plots for all samples

## Visualization

The pipeline generates two types of visualization:

1. **Coverage plots**: Shows the read depth across each segment with primer binding sites marked at the top
2. **Alignment visualizations**: Shows the alignment between original reference and metaconsensus

To view the visualizations, check the `results/plots/` directory.

## Troubleshooting

- Check log files in the sample directories for error messages
- Ensure input files follow the expected naming conventions
- Verify reference sequences are properly formatted FASTA files
- Make sure the conda environment is properly activated before running the pipeline

## License

[Include license information if applicable] 