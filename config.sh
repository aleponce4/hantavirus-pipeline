#!/bin/bash
#===============================================================================
# HANTAVIRUS SEQUENCING PIPELINE CONFIGURATION
#===============================================================================
# This file contains all configurable parameters for the pipeline.
# Edit these values to customize the pipeline behavior.
# The main pipeline will source this file to use these settings.

#-------------------------------------------------------------------------------
# GENERAL SETTINGS
#-------------------------------------------------------------------------------
# Whether to output detailed information during pipeline execution
VERBOSE=true

#-------------------------------------------------------------------------------
# QUALITY CONTROL & TRIMMING PARAMETERS
#-------------------------------------------------------------------------------
# Minimum quality score for read trimming (higher = more stringent)
TRIM_QUALITY=20

# Minimum read length to keep after trimming
MIN_READ_LENGTH=36

# Number of bases of overlap required for adapter detection
ADAPTER_STRINGENCY=3

#-------------------------------------------------------------------------------
# ALIGNMENT PARAMETERS
#-------------------------------------------------------------------------------
# Minimum seed length for BWA MEM (smaller = more sensitive but slower)
BWA_SEED_LENGTH=15

# Mismatch penalty for BWA MEM (higher = fewer mismatches allowed)
BWA_MISMATCH_PENALTY=3

#-------------------------------------------------------------------------------
# PRIMER HANDLING
#-------------------------------------------------------------------------------
# Minimum quality threshold for primer trimming
PRIMER_TRIM_QUALITY=50

# Minimum base quality score in primer trimming window
PRIMER_BASE_QUALITY=20

# Sliding window size for primer quality trimming
PRIMER_WINDOW_SIZE=4

#-------------------------------------------------------------------------------
# VARIANT CALLING PARAMETERS
#-------------------------------------------------------------------------------
# Minimum read depth required to call a variant
MIN_VARIANT_DEPTH=10

# Minimum variant frequency threshold (as a decimal, e.g., 0.03 = 3%)
MIN_VARIANT_FREQ=0.03

# Minimum quality score for a variant to be considered
MIN_VARIANT_QUALITY=20

#-------------------------------------------------------------------------------
# CONSENSUS GENERATION PARAMETERS
#-------------------------------------------------------------------------------
# Minimum coverage required to include a position in consensus
# Positions with coverage below this will be masked with 'N'
MIN_CONSENSUS_COVERAGE=10

# Character used to mask low-coverage regions in consensus
MASK_CHAR="N"

#-------------------------------------------------------------------------------
# SAMPLE CLASSIFICATION PARAMETERS
#-------------------------------------------------------------------------------
# Coverage threshold for negative sample detection
# Samples with average coverage below this are considered negative
NEGATIVE_SAMPLE_THRESHOLD=50

#-------------------------------------------------------------------------------
# PRIMER MISMATCH DETECTION PARAMETERS
#-------------------------------------------------------------------------------
# Whether to perform primer mismatch detection
DETECT_PRIMER_MISMATCHES=true

# Minimum quality score for filtered variants used in primer mismatch detection
PRIMER_MISMATCH_MIN_QUALITY=20

# Minimum frequency for filtered variants used in primer mismatch detection
PRIMER_MISMATCH_MIN_FREQ=0.03

#-------------------------------------------------------------------------------
# METACONSENSUS PARAMETERS
#-------------------------------------------------------------------------------
# Minimum sample count required to include a position in metaconsensus
MIN_SAMPLE_COUNT=2

#-------------------------------------------------------------------------------
# VARIANT ANNOTATION PARAMETERS
#-------------------------------------------------------------------------------
# Whether to perform enhanced variant annotation using SnpEff
ANNOTATE_VARIANTS=true

# Directory to install and run SnpEff (created automatically if needed)
SNPEFF_DIR="./tools/snpEff"

# Version of SnpEff to download
SNPEFF_VERSION="5.1"

# Java memory allocation for SnpEff
SNPEFF_MEMORY="2G"

# Reference IDs for segments (should match GenBank accession)
S_SEGMENT_REFID="JN232078"  # Update with your S segment reference ID
M_SEGMENT_REFID="OR184986"  # Update with your M segment reference ID
L_SEGMENT_REFID="JN232080"  # Update with your L segment reference ID 