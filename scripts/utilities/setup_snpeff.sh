#!/bin/bash
# Script to download and set up SnpEff for variant annotation

# Source the configuration file
source ./config.sh

# Store the original working directory
ORIGINAL_DIR=$(pwd)

# Check if SnpEff directory already exists
if [ -d "$SNPEFF_DIR" ] && [ -f "$SNPEFF_DIR/snpEff.jar" ]; then
    echo "SnpEff installation found at $SNPEFF_DIR"
else
    echo "Setting up SnpEff for variant annotation..."
    
    # Create the tools directory if it doesn't exist
    mkdir -p tools
    
    # Download SnpEff
    echo "Downloading SnpEff..."
    wget -q -O tools/snpEff_v${SNPEFF_VERSION}_core.zip https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
    
    # Extract SnpEff
    echo "Extracting SnpEff..."
    unzip -q -o tools/snpEff_v${SNPEFF_VERSION}_core.zip -d tools/
    
    # Rename directory if needed
    if [ ! -d "$SNPEFF_DIR" ]; then
        mv tools/snpEff "$SNPEFF_DIR"
    fi
    
    # Clean up zip file
    rm tools/snpEff_v${SNPEFF_VERSION}_core.zip
    
    echo "SnpEff installed successfully in $SNPEFF_DIR"
fi

# Create custom Hantavirus database for each segment
echo "Setting up custom Hantavirus database for SnpEff..."

# Ensure we're in the SnpEff directory for commands
cd "$SNPEFF_DIR"

# Setup directories for each segment
mkdir -p "data/hanta_S_segment" "data/hanta_M_segment" "data/hanta_L_segment"

# Create genomes file entries (if not already in config)
if ! grep -q "hanta_S_segment.genome" "snpEff.config"; then
    echo "Adding Hantavirus genomes to SnpEff configuration..."
    
    # Add configuration lines
    echo "# Hantavirus genomes" >> "snpEff.config"
    echo "hanta_S_segment.genome : Hantavirus S Segment" >> "snpEff.config"
    echo "hanta_M_segment.genome : Hantavirus M Segment" >> "snpEff.config"
    echo "hanta_L_segment.genome : Hantavirus L Segment" >> "snpEff.config"
fi

# Make sure we have reference genome and GFF files for each segment
for segment in "S" "M" "L"; do
    segment_dir="$ORIGINAL_DIR/data/references/${segment}_segment"
    ref_var="${segment}_SEGMENT_REFID"
    ref_id=${!ref_var}
    
    # Define the variable name we want to check
    segment_lower=$(echo "$segment" | tr '[:upper:]' '[:lower:]')
    genome_dir="data/hanta_${segment}_segment"
    
    # Skip if no reference ID is defined
    if [ -z "$ref_id" ]; then
        echo "Warning: No reference ID defined for ${segment} segment"
        continue
    fi
    
    # Copy and adapt reference FASTA
    if [ -f "$segment_dir/metaconsensus.fasta" ]; then
        # Use metaconsensus if available but change the header to match GFF ID
        echo "Using metaconsensus for ${segment} segment SnpEff database (with modified header)"
        # Create a temporary file with modified header
        sed "1s/^>metaconsensus/>$ref_id.1/" "$segment_dir/metaconsensus.fasta" > "$genome_dir/sequences.fa"
    elif [ -f "$segment_dir/${ref_id}.fasta" ]; then
        # Use reference if available and make sure header uses the correct ID
        echo "Using reference for ${segment} segment SnpEff database"
        # Create a temporary file with proper header
        sed "1s/^>.*/>$ref_id.1/" "$segment_dir/${ref_id}.fasta" > "$genome_dir/sequences.fa"
    else
        echo "Warning: No reference found for ${segment} segment"
        continue
    fi
    
    # Copy GFF file
    if [ -f "$segment_dir/${ref_id}.simple.gff" ]; then
        # Copy the reference GFF
        cp "$segment_dir/${ref_id}.simple.gff" "$genome_dir/genes.gff"
    elif [ -f "$segment_dir/metaconsensus.simple.gff" ]; then
        # Copy the metaconsensus GFF but update the seqid field if necessary
        echo "Using metaconsensus GFF with adjusted sequence ID"
        sed "s/^metaconsensus/$ref_id.1/" "$segment_dir/metaconsensus.simple.gff" > "$genome_dir/genes.gff"
    else
        echo "Warning: No GFF file found for ${segment} segment"
        continue
    fi
    
    # For L segment which has no data, skip build
    if [ "$segment" == "L" ]; then
        echo "Skipping database build for L segment (no data available)"
        continue
    fi

    # Try a simpler approach to build the database - just create a minimal predictor file
    echo "Building SnpEff database for ${segment} segment (simplified approach)..."
    
    # Create a simplified database structure
    mkdir -p "$genome_dir/regulation"
    touch "$genome_dir/regulation/regulation.gff"
    
    # Create a minimal predictor file directly
    # This tells SnpEff to just use the specified gene and exon regions without extensive validation
    PREDICTOR_DIR="$genome_dir"
    PREDICTOR_FILE="$PREDICTOR_DIR/snpEffectPredictor.bin"
    
    echo "Creating minimal database files for ${segment} segment..."
    
    # Force SnpEff to build the database with minimal checks
    java -Xmx${SNPEFF_MEMORY} -jar snpEff.jar build -noCheckCds -noCheckProtein -noLog -noDownload \
        -dataDir ./data -gff3 -v "hanta_${segment}_segment"
    
    # Verify if the predictor file was created
    if [ -f "$PREDICTOR_FILE" ]; then
        echo "Successfully created SnpEff database for ${segment} segment"
    else
        echo "WARNING: Failed to create SnpEff database for ${segment} segment"
    fi
done

# Go back to original directory
cd "$ORIGINAL_DIR"

echo "SnpEff setup completed." 