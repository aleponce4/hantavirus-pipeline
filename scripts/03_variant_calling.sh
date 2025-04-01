#!/bin/bash
# Variant calling script using LoFreq (parallel mode optimized for small viral segments)

sample=$1
threads=4  

if [ "$FIRST_PASS" = "true" ]; then
    RESULTS_DIR="results/first_pass"
else
    RESULTS_DIR="results/second_pass"
fi

for segment in L_segment M_segment S_segment; do
    segment_dir="data/references/$segment"
    if [ -d "$segment_dir" ] && [ "$(ls -A "$segment_dir")" ]; then

        if [ "$FIRST_PASS" = "true" ]; then
            reference=$(ls "$segment_dir"/*.fasta | grep -v "metaconsensus.fasta" | head -n 1)
        else
            metaconsensus="$segment_dir/metaconsensus.fasta"
            if [ ! -f "$metaconsensus" ]; then
                echo "Metaconsensus not found for $segment"
                exit 1
            fi
            reference="$metaconsensus"
        fi

        if [ ! -f "${reference}.fai" ]; then
            samtools faidx "$reference"
        fi

        bam_file="$RESULTS_DIR/$sample/alignment/$segment/${sample}.bam"
        out_dir="$RESULTS_DIR/$sample/variants/$segment"
        mkdir -p "$out_dir"

        rm -f "$out_dir/${sample}.vcf" "$out_dir/${sample}.vcf.gz"

        echo "Calling variants with LoFreq for $sample - $segment"
        lofreq call-parallel --pp-threads "$threads" --call-indels \
            -f "$reference" -o "$out_dir/${sample}.vcf" "$bam_file"

        bgzip -c "$out_dir/${sample}.vcf" > "$out_dir/${sample}.vcf.gz"
        tabix -p vcf "$out_dir/${sample}.vcf.gz"

        echo "Variants called and indexed for $sample - $segment"
    fi
done

echo "Variant calling completed for sample: $sample"
