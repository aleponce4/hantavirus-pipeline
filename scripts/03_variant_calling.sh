#!/bin/bash
# Full iVar pipeline with realignment to consensus for primer mismatch detection

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
        reference=$(ls "$segment_dir"/Reference_*.fasta | head -n 1)
        echo "Using reference for $segment: $reference"
        
        # Ensure reference is indexed
        if [ ! -f "${reference}.fai" ]; then
            samtools faidx "$reference"
            bwa index "$reference"
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
        
        # Step 2: Generate consensus genome with iVar
        echo "  Creating consensus genome..."
        samtools mpileup -aa -A -d 0 -Q 0 --reference "$reference" "$working_bam" | \
            ivar consensus -p "$out_dir/consensus" -m $MIN_CONSENSUS_COVERAGE -t 0.5 -q $MIN_VARIANT_QUALITY
        
        # Check if consensus was generated
        if [ ! -f "$out_dir/consensus.fa" ]; then
            echo "  ERROR: Failed to generate consensus for $sample - $segment"
            exit 1
        fi
        
        # Index the consensus for realignment
        samtools faidx "$out_dir/consensus.fa"
        bwa index "$out_dir/consensus.fa"
        
        # Step 3: REALIGN original reads to the consensus genome (FULL IVAR APPROACH)
        echo "  Realigning reads to consensus genome..."
        
        # Extract original FASTQ files from the initial BAM (before primer trimming)
        echo "  Extracting reads for realignment..."
        samtools sort -n "$bam_file" | \
            samtools fastq -@ $threads -1 "$out_dir/R1.fastq" -2 "$out_dir/R2.fastq" \
            -0 /dev/null -s "$out_dir/singleton.fastq"
        
        # Realign to consensus
        echo "  Mapping reads to consensus..."
        bwa mem -t $threads "$out_dir/consensus.fa" "$out_dir/R1.fastq" "$out_dir/R2.fastq" | \
            samtools view -bS - | \
            samtools sort -o "$out_dir/realigned_to_consensus.bam" -
        
        # Handle singleton reads if any
        if [ -s "$out_dir/singleton.fastq" ]; then
            bwa mem -t $threads "$out_dir/consensus.fa" "$out_dir/singleton.fastq" | \
                samtools view -bS - | \
                samtools sort -o "$out_dir/realigned_singleton.bam" -
            
            # Merge paired and singleton alignments
            samtools merge -f "$out_dir/realigned_merged.bam" "$out_dir/realigned_to_consensus.bam" "$out_dir/realigned_singleton.bam"
            mv "$out_dir/realigned_merged.bam" "$out_dir/realigned_to_consensus.bam"
            rm -f "$out_dir/realigned_singleton.bam"
        fi
        
        samtools index "$out_dir/realigned_to_consensus.bam"
        
        # Verify realignment worked
        read_count=$(samtools view -c "$out_dir/realigned_to_consensus.bam")
        echo "  Realigned BAM contains $read_count reads"
        
        # Step 4: Trim primers on realigned reads
        if [ -f "$primer_bed" ]; then
            echo "  Trimming primers on realigned reads..."
            ivar trim -i "$out_dir/realigned_to_consensus.bam" -b "$primer_bed" -p "$out_dir/realigned_trimmed" \
                -m $PRIMER_TRIM_QUALITY -q $PRIMER_BASE_QUALITY -s $PRIMER_WINDOW_SIZE -e
            
            samtools sort -o "$out_dir/realigned_trimmed.sorted.bam" "$out_dir/realigned_trimmed.bam"
            samtools index "$out_dir/realigned_trimmed.sorted.bam"
            final_bam="$out_dir/realigned_trimmed.sorted.bam"
        else
            final_bam="$out_dir/realigned_to_consensus.bam"
        fi
        
        # Step 5: Detect primer mismatches using realigned data (FULL IVAR APPROACH)
        if [ "$DETECT_PRIMER_MISMATCHES" = "true" ] && [ -f "$primer_bed" ] && [ -f "$primer_pairs" ]; then
            echo "  Checking for primer mismatches using realigned data..."
            
            # Call variants against consensus using realigned reads
            samtools mpileup -aa -A -d 0 -Q 0 --reference "$out_dir/consensus.fa" "$final_bam" | \
                ivar variants -p "$out_dir/primer_check" -r "$out_dir/consensus.fa" \
                -m $MIN_VARIANT_DEPTH -t $MIN_VARIANT_FREQ $gff_arg
            
            # Identify problematic primers
            if [ -f "$out_dir/primer_check.tsv" ]; then
                ivar getmasked -i "$out_dir/primer_check.tsv" -b "$primer_bed" -f "$primer_pairs" -p "$out_dir/masked"
                
                if [ -f "$out_dir/masked.txt" ] && [ -s "$out_dir/masked.txt" ]; then
                    masked_count=$(wc -w < "$out_dir/masked.txt")
                    echo "  Found $masked_count problematic primers"
                    
                    # Remove reads from problematic amplicons
                    echo "  Removing reads from problematic amplicons..."
                    ivar removereads -i "$final_bam" -p "$out_dir/cleaned" -t "$out_dir/masked.txt" -b "$primer_bed"
                    
                    samtools sort -o "$out_dir/cleaned.sorted.bam" "$out_dir/cleaned.bam"
                    samtools index "$out_dir/cleaned.sorted.bam"
                    final_bam="$out_dir/cleaned.sorted.bam"
                    echo "  Using cleaned BAM for final analysis"
                else
                    echo "  No problematic primers detected"
                fi
            fi
        else
            echo "  Primer mismatch detection disabled or files missing"
        fi
        
        # Step 6: Final variant calling against ORIGINAL reference using cleaned realigned data
        echo "  Final variant calling against original reference..."
        
        # For final variants, we need to realign the cleaned reads back to original reference
        # Extract reads from cleaned BAM
        samtools sort -n "$final_bam" | \
            samtools fastq -@ $threads -1 "$out_dir/final_R1.fastq" -2 "$out_dir/final_R2.fastq" \
            -0 /dev/null -s "$out_dir/final_singleton.fastq"
        
        # Realign to original reference
        bwa mem -t $threads "$reference" "$out_dir/final_R1.fastq" "$out_dir/final_R2.fastq" | \
            samtools view -bS - | \
            samtools sort -o "$out_dir/final_alignment.bam" -
        
        if [ -s "$out_dir/final_singleton.fastq" ]; then
            bwa mem -t $threads "$reference" "$out_dir/final_singleton.fastq" | \
                samtools view -bS - | \
                samtools sort -o "$out_dir/final_singleton.bam" -
            
            samtools merge -f "$out_dir/final_merged.bam" "$out_dir/final_alignment.bam" "$out_dir/final_singleton.bam"
            mv "$out_dir/final_merged.bam" "$out_dir/final_alignment.bam"
            rm -f "$out_dir/final_singleton.bam"
        fi
        
        samtools index "$out_dir/final_alignment.bam"
        
        # Call variants against original reference
        samtools mpileup -aa -A -d 0 -Q 0 --reference "$reference" "$out_dir/final_alignment.bam" | \
            ivar variants -p "$out_dir/${sample}" -r "$reference" \
            -m $MIN_VARIANT_DEPTH -t $MIN_VARIANT_FREQ $gff_arg
        
        # Convert to VCF for compatibility
        echo "  Converting to VCF format..."
        python scripts/utilities/ivar_variants_to_vcf.py "$out_dir/${sample}.tsv" "$reference" > "$out_dir/${sample}.vcf"
        bgzip -c "$out_dir/${sample}.vcf" > "$out_dir/${sample}.vcf.gz"
        tabix -p vcf "$out_dir/${sample}.vcf.gz"
        
        # Report results
        variant_count=$(grep -v "^REGION" "$out_dir/${sample}.tsv" | wc -l)
        echo "  Found $variant_count variants"
        
        if [ -f "$gff_file" ]; then
            nonsyn_count=$(grep -v "Synonymous" "$out_dir/${sample}.tsv" | grep -v "REGION" | grep -v "REF" | wc -l)
            echo "  Including $nonsyn_count non-synonymous mutations"
        fi
        
        # Clean up intermediate files
        rm -f "$out_dir/trimmed.bam" "$out_dir/realigned_trimmed.bam" "$out_dir/cleaned.bam"
        rm -f "$out_dir/R1.fastq" "$out_dir/R2.fastq" "$out_dir/singleton.fastq"
        rm -f "$out_dir/final_R1.fastq" "$out_dir/final_R2.fastq" "$out_dir/final_singleton.fastq"
        
        echo "  Completed $sample - $segment"
    fi
done

echo "Full iVar pipeline completed for sample: $sample"
