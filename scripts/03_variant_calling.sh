#!/bin/bash
# Simplified variant calling script using iVar with primer mismatch detection

sample=$1
threads=${LOFREQ_THREADS:-4}

# Source the configuration file
source ./config.sh

if [ "$FIRST_PASS" = "true" ]; then
    RESULTS_DIR="results/first_pass"
    echo "Processing $sample (FIRST PASS)"
else
    RESULTS_DIR="results/second_pass"
    echo "Processing $sample (SECOND PASS)"
fi

# Create primer BED files if they don't exist
if [ ! -d "data/primers/bed_files" ] || [ ! "$(ls -A data/primers/bed_files)" ]; then
    echo "Generating primer BED files..."
    bash scripts/utilities/generate_primer_bed.sh
fi

# Create primer pair files if they don't exist
if [ ! -d "data/primers/pairings" ] || [ ! "$(ls -A data/primers/pairings)" ]; then
    echo "Generating primer pair files..."
    python scripts/utilities/generate_primer_pairs.py data/primers/Primers.csv data/primers/pairings
fi

# Generate simplified GFF files for iVar amino acid annotation
if [ ! -f "data/references/S_segment/Reference_S.simple.gff" ] || [ ! -f "data/references/M_segment/Reference_M.simple.gff" ]; then
    echo "Generating simplified GFF files for amino acid annotation..."
    bash scripts/utilities/create_simple_gff.sh
fi

for segment in M_segment S_segment; do  # Skip L_segment - no data
    segment_dir="data/references/$segment"
    
    if [ -d "$segment_dir" ] && [ "$(ls -A "$segment_dir")" ]; then
        # Reference selection: use original reference for both passes
        # (metaconsensus complexity removed)
        reference=$(ls "$segment_dir"/Reference_*.fasta | head -n 1)
        echo "Using reference for $segment: $reference"
        
        # Ensure reference is indexed
        if [ ! -f "${reference}.fai" ]; then
            samtools faidx "$reference"
        fi
        
        # Set up paths
        bam_file="$RESULTS_DIR/$sample/alignment/$segment/${sample}.bam"
        out_dir="$RESULTS_DIR/$sample/variants/$segment"
        mkdir -p "$out_dir"
        
        # Set up primer files
        primer_bed="data/primers/bed_files/${segment}_primers.bed"
        segment_short=$(echo "$segment" | cut -d'_' -f1)
        primer_pairs="data/primers/pairings/${segment_short}_segment_primer_pairs.tsv"
        
        # Set up GFF for amino acid annotation
        gff_file="data/references/$segment/Reference_${segment_short}.simple.gff"
        if [ -f "$gff_file" ]; then
            gff_arg="-g $gff_file"
            echo "Using GFF file for amino acid annotations: $gff_file"
        else
            gff_arg=""
            echo "No GFF file found - amino acid annotations will not be available"
        fi
        
        echo "Processing $sample - $segment"
        
        # Step 1: Trim primers with iVar
        if [ -f "$primer_bed" ]; then
            echo "  Trimming primers with iVar..."
            ivar trim -i "$bam_file" -b "$primer_bed" -p "$out_dir/trimmed" \
                -m $PRIMER_TRIM_QUALITY -q $PRIMER_BASE_QUALITY -s $PRIMER_WINDOW_SIZE -e
            
            samtools sort -o "$out_dir/trimmed.sorted.bam" "$out_dir/trimmed.bam"
            samtools index "$out_dir/trimmed.sorted.bam"
            working_bam="$out_dir/trimmed.sorted.bam"
            echo "  Using primer-trimmed BAM"
        else
            echo "  No primer BED file found - proceeding without primer trimming"
            working_bam="$bam_file"
        fi
        
        # Step 2: Generate draft consensus with iVar
        echo "  Creating draft consensus..."
        samtools mpileup -aa -A -d 0 -Q 0 --reference "$reference" "$working_bam" | \
            ivar consensus -p "$out_dir/draft_consensus" -m $MIN_CONSENSUS_COVERAGE -t 0.5 -q $MIN_VARIANT_QUALITY
        
        # Step 3: Detect primer mismatches (if enabled and files exist)
        if [ "$DETECT_PRIMER_MISMATCHES" = "true" ] && [ -f "$primer_bed" ] && [ -f "$primer_pairs" ]; then
            echo "  Checking for primer mismatches..."
            
            # Call variants against draft consensus to check primer binding
            samtools mpileup -aa -A -d 0 -Q 0 --reference "$out_dir/draft_consensus.fa" "$working_bam" | \
                ivar variants -p "$out_dir/primer_check" -r "$out_dir/draft_consensus.fa" \
                -m $MIN_VARIANT_DEPTH -t $MIN_VARIANT_FREQ $gff_arg
            
            # Identify problematic primers
            if [ -f "$out_dir/primer_check.tsv" ]; then
                ivar getmasked -i "$out_dir/primer_check.tsv" -b "$primer_bed" -f "$primer_pairs" -p "$out_dir/masked"
                
                if [ -f "$out_dir/masked.txt" ] && [ -s "$out_dir/masked.txt" ]; then
                    masked_count=$(wc -w < "$out_dir/masked.txt")
                    echo "  Found $masked_count problematic primers"
                    
                    # Remove reads from problematic amplicons
                    echo "  Removing reads from problematic amplicons..."
                    ivar removereads -i "$working_bam" -p "$out_dir/cleaned" -t "$out_dir/masked.txt" -b "$primer_bed"
                    
                    samtools sort -o "$out_dir/cleaned.sorted.bam" "$out_dir/cleaned.bam"
                    samtools index "$out_dir/cleaned.sorted.bam"
                    working_bam="$out_dir/cleaned.sorted.bam"
                    echo "  Using cleaned BAM for final analysis"
                else
                    echo "  No problematic primers detected"
                fi
            fi
        else
            echo "  Primer mismatch detection disabled or files missing"
        fi
        
        # Step 4: Final variant calling with iVar
        echo "  Calling variants..."
        samtools mpileup -aa -A -d 0 -Q 0 --reference "$reference" "$working_bam" | \
            ivar variants -p "$out_dir/${sample}" -r "$reference" \
            -m $MIN_VARIANT_DEPTH -t $MIN_VARIANT_FREQ $gff_arg
        
        # Step 5: Convert to VCF for compatibility
        echo "  Converting to VCF format..."
        python scripts/utilities/ivar_variants_to_vcf.py "$out_dir/${sample}.tsv" "$reference" > "$out_dir/${sample}.vcf"
        bgzip -c "$out_dir/${sample}.vcf" > "$out_dir/${sample}.vcf.gz"
        tabix -p vcf "$out_dir/${sample}.vcf.gz"
        
        # Report results
        variant_count=$(grep -v "^#" "$out_dir/${sample}.tsv" | wc -l)
        echo "  Found $variant_count variants"
        
        if [ -f "$gff_file" ]; then
            nonsyn_count=$(grep -v "Synonymous" "$out_dir/${sample}.tsv" | grep -v "REGION" | grep -v "REF" | wc -l)
            echo "  Including $nonsyn_count non-synonymous mutations"
        fi
        
        # Clean up intermediate files
        rm -f "$out_dir/trimmed.bam" "$out_dir/cleaned.bam" "$out_dir/draft_consensus.qual.txt"
        
        echo "  Completed $sample - $segment"
    fi
done

echo "Variant calling completed for sample: $sample"
