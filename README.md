# Hantavirus Pipeline

A pipeline for processing and analyzing Hantavirus NGS data, including quality control, mapping, variant calling, consensus generation, and visualization.

## Overview

This pipeline is designed for processing Illumina paired-end sequencing data of Hantavirus genomes. It features a two-pass approach with automatic detection of negative samples based on coverage depth. The pipeline handles:

- Quality trimming and adapter removal
- Primer trimming
- Reference-based mapping
- Automatic negative sample detection
- Metaconsensus generation
- Variant calling with improved coordinate consistency
- Consensus sequence creation
- Coverage visualization
- Primer evaluation and optimization

## Installation

### Prerequisites

- Conda or Miniconda (download from [https://docs.conda.io/projects/miniconda/en/latest/](https://docs.conda.io/projects/miniconda/en/latest/))
- Git (optional, for cloning the repository)

### Setup

1. Clone or download the repository:
   ```bash
   git clone [https://github.com/aleponce4/hantavirus-pipeline]
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
   - GFF3 annotation files should be named to match the reference FASTA (e.g., `JN232078.fasta` and `JN232078.gff3`)

3. Place your primer information in `data/primers/Primers.csv` (if using primer trimming)

4. Primer File Format (data/primers/Primers.csv)
   The pipeline expects a CSV file with the following columns:
   - Column 1: Primer Name - Using the following naming convention:
     * SF* or SR* for S segment primers (forward and reverse)
     * MF* or MR* for M segment primers (forward and reverse)
     * LF* or LR* for L segment primers (forward and reverse)
   - Column 2: Region - Genomic position in format like "1-32F" or "452-480R"
   - Column 3: Sequence - The primer sequence (spaces will be removed)
   - Additional columns may be present but are not required for basic functionality

   Example header line: 
   ```
   F Name,Region,Sequence,Tm,GC%,Amplicon
   ```

   Example data line:
   ```
   SF1,1-32F,TAGTAGTAGACTCCTTGAGAAGCTACT,58.2,40.7,S1
   ```

### Running the Pipeline

1. Activate the conda environment:
   ```bash
   conda activate De_Novo_pipeline
   ```

2. Run the pipeline:
   ```bash
   bash run_pipeline.sh
   ```

3. Alternatively, to run multiple samples sequentially:
   ```bash
   for sample in sample1 sample2 sample3; do
       bash run_pipeline.sh
   done
   ```

### HPC Usage

The pipeline includes optimizations for high-performance computing (HPC) environments:

1. **Automatic Resource Detection**:
   - The pipeline automatically detects systems with >24 cores and enables HPC mode
   - Thread allocation is optimized for different tools based on their scaling characteristics

2. **Sample-Level Parallelism**:
   - On HPC systems, multiple samples will be processed concurrently
   - Default configuration processes 3 samples simultaneously on a 48-core system
   - Each sample uses optimized thread counts for specific tools:
     * BWA: 16 threads per sample
     * Samtools: 8 threads per sample 
     * LoFreq: 8 threads per sample
     * Trim Galore: 8 threads per sample
     * MAFFT: All available threads (during metaconsensus generation)

3. **Manual Configuration**:
   - Edit the pipeline.sh file to adjust MAX_CONCURRENT_SAMPLES and thread counts
   - For systems with more RAM, you can increase concurrent samples
   - For memory-limited systems, reduce threads per tool

4. **Job Management**:
   - The pipeline manages subprocesses automatically
   - Limits concurrency to prevent system overload
   - Waits for jobs to complete before moving to the next pipeline phase
   - Creates job tracking markers in the .jobs/ directory

5. **HPC Job Submission Example**:
   ```bash
   # Example SLURM job submission script
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

### Pipeline Steps

The pipeline executes the following steps:

1. **First Pass Processing**:
   - Quality trimming, adapter removal, and primer trimming
   - Alignment to original reference sequences
   - Variant calling
   - Consensus generation

2. **Negative Sample Detection**:
   - Calculates average coverage for each segment
   - Samples with average coverage below 50x are classified as negative
   - Negative samples skip second pass processing

3. **Reference Refinement**:
   - Creates a metaconsensus sequence from all positive samples

4. **Second Pass Processing** (positive samples only):
   - Re-alignment to metaconsensus reference (for better read mapping)
   - Extracting reads from the metaconsensus alignment
   - Re-mapping extracted reads to the original reference
   - Variant calling against the original reference (for consistent coordinates)
   - Final consensus generation
   - Comparison between first and second pass results

5. **Visualization and Analysis**:
   - Coverage plots for positive and negative samples
   - Alignment visualizations
   - Consensus comparison metrics
   
6. **Primer Evaluation**:
   - Evaluates primer quality against consensus sequences
   - Generates improvement suggestions for problematic primers
   - Creates primer maps and visualizations

### Output Directory Structure

Results are organized in the `results/` directory with the following structure:

- `results/first_pass/` - First pass processing results for all samples
- `results/second_pass/` - Second pass processing results (positive samples only)
- `results/negative_samples/` - Data for samples with low coverage
- `results/trimmed/` - Trimmed reads (used by both passes)
- `results/plots/` - Coverage and alignment plots
  - `results/plots/negative_samples/` - Plots for negative samples
- `results/primer_evaluation/` - Primer evaluation results and suggested improvements

## Variant Calling Improvements

The pipeline implements a sophisticated two-pass approach to variant calling:

1. In the first pass, reads are mapped directly to the original reference
2. In the second pass:
   - Reads are initially mapped to the metaconsensus (better mapping due to closer sequence similarity)
   - Mapped reads are extracted and remapped to the original reference
   - Variants are called against the original reference
   
This approach solves several challenges:
- Improves read mapping by using the metaconsensus (which better represents the virus population)
- Maintains coordinate consistency by calling variants against the original reference
- Prevents false positive variants caused by coordinate mismatches
- Enables proper amino acid annotation using consistent GFF coordinates

## Visualization

The pipeline generates multiple types of visualizations:

1. **Coverage plots**: Shows the read depth across each segment with primer binding sites marked
2. **Alignment visualizations**: Shows the alignment between original reference and metaconsensus
3. **Primer maps**: Visualizations of primer binding locations and quality metrics

## Primer Evaluation

The pipeline includes a primer evaluation module that:

1. Analyzes primer quality metrics (Tm, GC content, hairpins, etc.)
2. Evaluates primer binding efficiency to consensus sequences
3. Identifies problematic primers
4. Generates improved primer sequences when possible
5. Creates summary reports and visualizations

## Troubleshooting

- Check log files in the `logs/` directory for error messages
- Ensure input files follow the expected naming conventions
- Verify reference sequences are properly formatted FASTA files
- Make sure the conda environment is properly activated before running the pipeline
- For negative samples, check `results/negative_samples/negative_samples.txt` for the list of samples that were classified as negative
- When using custom references, ensure that GFF files match the FASTA files for proper annotation

## License

This pipeline is provided for research use. Please cite [Citation information] when using this pipeline in your research.