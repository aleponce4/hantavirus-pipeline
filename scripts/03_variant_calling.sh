#!/bin/bash
# Variant calling script using iVar with primer masking and amino acid annotations

sample=$1
threads=${LOFREQ_THREADS:-4}  # Use environment variable if set, otherwise default to 4

# Source the configuration file
source ./config.sh

if [ "$FIRST_PASS" = "true" ]; then
    RESULTS_DIR="results/first_pass"
else
    RESULTS_DIR="results/second_pass"
fi

# Create primer BED files if they don't exist yet
if [ ! -d "data/primers/bed_files" ] || [ ! "$(ls -A data/primers/bed_files)" ]; then
    echo "Generating primer BED files..."
    bash scripts/utilities/generate_primer_bed.sh
fi

# Create primer pair files if they don't exist yet
if [ ! -d "data/primers/pairings" ] || [ ! "$(ls -A data/primers/pairings)" ]; then
    echo "Generating primer pair files..."
    python scripts/utilities/generate_primer_pairs.py data/primers/Primers.csv data/primers/pairings
fi

# Generate simplified GFF files for annotation
if [ ! -f "data/references/S_segment/JN232078.simple.gff" ] || [ ! -f "data/references/M_segment/OR184986.simple.gff" ]; then
    echo "Generating simplified GFF files for amino acid annotation..."
    bash scripts/utilities/create_simple_gff.sh
fi

# If we're doing the second pass, we need GFF files for metaconsensus
if [ "$FIRST_PASS" != "true" ]; then
    # Generate metaconsensus GFF files if they don't exist
    if [ ! -f "data/references/S_segment/metaconsensus.simple.gff" ] || [ ! -f "data/references/M_segment/metaconsensus.simple.gff" ]; then
        echo "Generating metaconsensus GFF files for second pass amino acid annotation..."
        bash scripts/utilities/create_metaconsensus_gff.sh
    fi
fi

# If we want to use SnpEff for annotation, ensure it's set up
if [ "$ANNOTATE_VARIANTS" = "true" ] && [ ! -f "$SNPEFF_DIR/snpEff.jar" ]; then
    echo "Setting up SnpEff for variant annotation..."
    bash scripts/utilities/setup_snpeff.sh
fi

for segment in L_segment M_segment S_segment; do
    segment_dir="data/references/$segment"
    if [ -d "$segment_dir" ] && [ "$(ls -A "$segment_dir")" ]; then

        # For the first pass, use the original reference for alignment and variant calling
        if [ "$FIRST_PASS" = "true" ]; then
            reference=$(ls "$segment_dir"/*.fasta | grep -v "metaconsensus.fasta" | head -n 1)
            variant_reference=$reference
            gff_reference=$(basename "$reference" .fasta)
        else
            # For second pass, use metaconsensus for alignment but original reference for variant calling
            metaconsensus="$segment_dir/metaconsensus.fasta"
            if [ ! -f "$metaconsensus" ]; then
                echo "Metaconsensus not found for $segment"
                exit 1
            fi
            
            # The alignment was done against metaconsensus
            reference="$metaconsensus"
            
            # But for variant calling against the original reference
            variant_reference=$(ls "$segment_dir"/*.fasta | grep -v "metaconsensus.fasta" | head -n 1)
            gff_reference=$(basename "$variant_reference" .fasta)
        fi

        if [ ! -f "${reference}.fai" ]; then
            samtools faidx "$reference"
        fi
        
        if [ ! -f "${variant_reference}.fai" ]; then
            samtools faidx "$variant_reference"
        fi

        # Check for GFF file for amino acid annotations - use original reference GFF for consistent coordinates
        simple_gff="${segment_dir}/${gff_reference}.simple.gff"
        if [ -f "$simple_gff" ]; then
            echo "Found simplified GFF file for amino acid annotations: $simple_gff"
            gff_arg="-g $simple_gff"
        else
            echo "No simplified GFF file found for $segment. Amino acid changes will not be annotated."
            gff_arg=""
        fi

        bam_file="$RESULTS_DIR/$sample/alignment/$segment/${sample}.bam"
        out_dir="$RESULTS_DIR/$sample/variants/$segment"
        mkdir -p "$out_dir"

        # Define the bed file for primer masking
        primer_bed="data/primers/bed_files/${segment}_primers.bed"
        
        # Define the primer pair file
        segment_short=$(echo "$segment" | cut -d'_' -f1)
        primer_pairs="data/primers/pairings/${segment_short}_segment_primer_pairs.tsv"

        echo "Calling variants for $sample - $segment"
        
        # Step 1: Trim primers using iVar if primer BED file exists
        if [ -f "$primer_bed" ]; then
            echo "  Masking primers using iVar trim..."
            # Create a temporary directory for intermediate files
            tmp_dir=$(mktemp -d)
            
            # Use iVar to trim primers
            ivar trim -i "$bam_file" -b "$primer_bed" -p "$tmp_dir/trimmed" \
                -m $PRIMER_TRIM_QUALITY -q $PRIMER_BASE_QUALITY -s $PRIMER_WINDOW_SIZE -e
            
            # Sort and index the trimmed BAM
            samtools sort -o "$tmp_dir/trimmed.sorted.bam" "$tmp_dir/trimmed.bam"
            samtools index "$tmp_dir/trimmed.sorted.bam"
            
            # Use trimmed BAM for initial variant calling
            input_bam="$tmp_dir/trimmed.sorted.bam"
            echo "  Using primer-trimmed BAM for variant calling"
        else
            echo "  No primer BED file found for $segment, proceeding without primer masking"
            input_bam="$bam_file"
        fi
        
        # For second pass, we need to convert alignment from metaconsensus to original reference
        if [ "$FIRST_PASS" != "true" ] && [ "$reference" != "$variant_reference" ]; then
            echo "  Second pass: Converting alignment from metaconsensus to original reference..."
            
            # Create a directory for the remapped BAM
            remap_dir="$tmp_dir/remap"
            mkdir -p "$remap_dir"
            
            # Extract reads from the BAM file
            echo "  Extracting reads from metaconsensus-aligned BAM..."
            # Modified command to ensure we capture all reads properly
            samtools sort -n "$input_bam" | \
                samtools fastq -@ $threads -1 "$remap_dir/R1.fastq" -2 "$remap_dir/R2.fastq" \
                -0 /dev/null -s "$remap_dir/singleton.fastq"
            
            # Check if reads were extracted properly
            r1_count=$(grep -c "^@" "$remap_dir/R1.fastq" || echo 0)
            r2_count=$(grep -c "^@" "$remap_dir/R2.fastq" || echo 0)
            echo "  Extracted $r1_count read pairs for remapping"
            
            # Map reads to the original reference
            echo "  Remapping reads to original reference..."
            bwa mem -t $threads "$variant_reference" "$remap_dir/R1.fastq" "$remap_dir/R2.fastq" | \
                samtools view -bS - | \
                samtools sort -o "$remap_dir/remapped.bam" -
            
            # Also map singleton reads if any exist
            if [ -s "$remap_dir/singleton.fastq" ]; then
                echo "  Remapping singleton reads..."
                bwa mem -t $threads "$variant_reference" "$remap_dir/singleton.fastq" | \
                    samtools view -bS - | \
                    samtools sort -o "$remap_dir/remapped_singleton.bam" -
                
                # Merge the paired and singleton BAMs
                samtools merge -f "$remap_dir/merged.bam" "$remap_dir/remapped.bam" "$remap_dir/remapped_singleton.bam"
                samtools index "$remap_dir/merged.bam"
                input_bam="$remap_dir/merged.bam"
            else
                # Index the remapped BAM
                samtools index "$remap_dir/remapped.bam"
                input_bam="$remap_dir/remapped.bam"
            fi
            
            # Verify the remapped BAM has reads
            read_count=$(samtools view -c "$input_bam")
            echo "  Remapped BAM contains $read_count reads"
            
            if [ "$read_count" -lt 100 ]; then
                echo "  WARNING: Very few reads in remapped BAM. Using original primer-trimmed BAM for variant calling."
                input_bam="$tmp_dir/trimmed.sorted.bam"
            else
                echo "  Using remapped BAM aligned to original reference for variant calling"
            fi
        fi
        
        # Step 2: Generate pileup against the variant reference
        echo "  Creating pileup..."
        samtools mpileup -aa -A -d 0 -Q 0 --reference "$variant_reference" "$input_bam" > "$out_dir/${sample}.pileup"
        
        # Step 3: Call variants with iVar, including amino acid annotations if GFF file exists
        echo "  Calling variants with iVar..."
        cat "$out_dir/${sample}.pileup" | ivar variants -p "$out_dir/${sample}" \
            -r "$variant_reference" -m $MIN_VARIANT_DEPTH -t $MIN_VARIANT_FREQ $gff_arg
        
        # Step 4: Filter variants for primer mismatch detection
        echo "  Filtering variants for primer mismatch detection..."
        bash scripts/utilities/filter_variants.sh "$sample" "$segment" "$RESULTS_DIR" $MIN_VARIANT_QUALITY $MIN_VARIANT_FREQ
        
        # Step 5: Detect primer mismatches if primer files exist and enabled in config
        if [ "$DETECT_PRIMER_MISMATCHES" = "true" ] && [ -f "$primer_bed" ] && [ -f "$primer_pairs" ] && [ -f "$input_bam" ]; then
            echo "  Detecting primer mismatches using iVar getmasked..."
            ivar getmasked -i "$out_dir/${sample}.filtered.tsv" -b "$primer_bed" -f "$primer_pairs" -p "$out_dir/${sample}.masked"
            
            if [ -f "$out_dir/${sample}.masked.txt" ]; then
                masked_count=$(wc -w < "$out_dir/${sample}.masked.txt")
                if [ "$masked_count" -gt 0 ]; then
                    echo "  Found $masked_count primers with potential mismatches"
                    
                    # Step 6: Remove reads from mismatched primers
                    echo "  Removing reads from mismatched primers..."
                    ivar removereads -i "$input_bam" -p "$out_dir/${sample}.cleaned" -t "$out_dir/${sample}.masked.txt" -b "$primer_bed"
                    
                    # Use the cleaned BAM for final variant calling
                    samtools sort -o "$out_dir/${sample}.cleaned.sorted.bam" "$out_dir/${sample}.cleaned.bam"
                    samtools index "$out_dir/${sample}.cleaned.sorted.bam"
                    
                    # Re-generate pileup with cleaned BAM
                    echo "  Re-calling variants with cleaned BAM file..."
                    samtools mpileup -aa -A -d 0 -Q 0 --reference "$variant_reference" "$out_dir/${sample}.cleaned.sorted.bam" > "$out_dir/${sample}.cleaned.pileup"
                    
                    # Re-call variants with cleaned BAM
                    cat "$out_dir/${sample}.cleaned.pileup" | ivar variants -p "$out_dir/${sample}.final" \
                        -r "$variant_reference" -m $MIN_VARIANT_DEPTH -t $MIN_VARIANT_FREQ $gff_arg
                    
                    # Use final variants as the main output
                    mv "$out_dir/${sample}.final.tsv" "$out_dir/${sample}.tsv"
                    
                    echo "  Variant calling completed with primer mismatch correction"
                else
                    echo "  No primer mismatches detected"
                fi
            else
                echo "  Warning: Primer mismatch detection failed, using standard variants"
            fi
        else
            if [ "$DETECT_PRIMER_MISMATCHES" != "true" ]; then
                echo "  Primer mismatch detection disabled in configuration"
            else
                echo "  Skipping primer mismatch detection (missing required files)"
            fi
        fi
        
        # Step 7: Convert iVar TSV output to VCF format
        echo "  Converting iVar output to VCF..."
        python scripts/utilities/ivar_variants_to_vcf.py "$out_dir/${sample}.tsv" "$variant_reference" > "$out_dir/${sample}.vcf"
        
        # Step 8: Compress and index VCF for downstream use
        bgzip -c "$out_dir/${sample}.vcf" > "$out_dir/${sample}.vcf.gz"
        tabix -p vcf "$out_dir/${sample}.vcf.gz"
        
        # Step 9: Enhance variants with SnpEff annotation if enabled
        if [ "$ANNOTATE_VARIANTS" = "true" ]; then
            echo "  Enhancing variant annotations with SnpEff..."
            bash scripts/utilities/annotate_variants.sh "$sample" "$segment" "$RESULTS_DIR"
        fi
        
        # Clean up temporary files if they exist
        if [ -d "$tmp_dir" ]; then
            rm -rf "$tmp_dir"
        fi

        echo "Variants called and indexed for $sample - $segment"
        
        # Report non-synonymous variants if GFF file was used
        if [ -f "$simple_gff" ]; then
            nonsyn_count=$(grep -v "Synonymous" "$out_dir/${sample}.tsv" | grep -v "NA" | grep -v "REGION" | wc -l)
            if [ $nonsyn_count -gt 0 ]; then
                echo "  Found $nonsyn_count non-synonymous mutations"
            fi
        fi
    fi
done

echo "Variant calling completed for sample: $sample"
