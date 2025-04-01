# Primer Evaluation Module

This module evaluates primers used in the Hantavirus pipeline according to established quality criteria.

## Features

- Calculates GC content and melting temperature
- Predicts secondary structure formation (hairpins)
- Assesses specificity through BLAST E-values
- Evaluates primer length and composition
- Generates visual maps of primer positions
- Produces comprehensive quality reports

## Usage

### As part of the pipeline

Primer evaluation now runs automatically after the main pipeline completes:

```bash
./run_pipeline.sh
```

### As a standalone tool

To evaluate primers independently:

```bash
bash ./scripts/07_evaluate_primers.sh
```

Or directly using the Python script:

```bash
python scripts/primer_evaluation/evaluate_primers.py \
    --primers data/primers/Primers.csv \
    --segment S_segment \
    --blast_db data/blast_db/S_segment_blast_db \
    --consensus data/references/S_segment/metaconsensus.fasta \
    --output results/primer_evaluation
```

## Quality Criteria

Primers are evaluated against the following criteria:

1. **Primer Length**: 18-35 bases
2. **GC Content**: 40-60%
3. **Melting Temperature (Tm)**: 64-68°C
4. **Secondary Structure**: Gibbs Free Energy > -3 kcal/mol
5. **Specificity**: BLAST E-value ≤ 10^-4
6. **Base Composition**: No long runs of the same base (e.g., AAAAA) or repeats (e.g., ATATAT)

## Output

Results are stored in the `results/primer_evaluation` directory:

- `S_segment_primer_evaluation.csv`: Detailed analysis of each primer
- `S_segment_primer_map.png`: Visual representation of primer positions

## Implementation Details

The primer evaluation uses:

- **Melting Temperature Calculation**: Nearest-neighbor thermodynamics (Allawi & SantaLucia, 1997)
- **Secondary Structure Prediction**: seqfold package for ΔG calculations
- **Specificity Analysis**: BLAST search against segment-specific consensus sequences

## Adding New Primers

To evaluate new primers:

1. Add them to the `data/primers/Primers.csv` file following the existing format
2. Run the primer evaluation module

## Troubleshooting

- Ensure the conda environment is properly activated
- Check that BLAST databases are correctly created
- Verify the primer CSV file format matches the expected format 