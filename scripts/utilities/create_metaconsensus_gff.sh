#!/bin/bash

# Script to create GFF files for metaconsensus references

# For S segment
S_META="data/references/S_segment/metaconsensus.fasta"
S_META_GFF="data/references/S_segment/metaconsensus.simple.gff"

# For M segment
M_META="data/references/M_segment/metaconsensus.fasta"
M_META_GFF="data/references/M_segment/metaconsensus.simple.gff"

# Get the existing ORF coordinates from original GFF files
S_ORIGINAL_GFF="data/references/S_segment/JN232078.simple.gff"
M_ORIGINAL_GFF="data/references/M_segment/OR184986.simple.gff"

# S segment: Extract CDS information and create metaconsensus GFF
if [ -f "$S_META" ] && [ -f "$S_ORIGINAL_GFF" ]; then
    # Get the sequence ID from metaconsensus FASTA
    s_meta_id=$(grep ">" "$S_META" | head -1 | sed 's/>//' | cut -d ' ' -f1)
    
    if [ -z "$s_meta_id" ]; then
        echo "Error: Could not extract sequence ID from $S_META"
    else
        echo "Creating GFF for S segment metaconsensus (ID: $s_meta_id)..."
        
        # Extract the start and end positions from the original GFF
        s_start=$(grep "CDS" "$S_ORIGINAL_GFF" | awk '{print $4}')
        s_end=$(grep "CDS" "$S_ORIGINAL_GFF" | awk '{print $5}')
        
        if [ -n "$s_start" ] && [ -n "$s_end" ]; then
            # Create a GFF file for the metaconsensus
            cat > "$S_META_GFF" << EOF
##gff-version 3
$s_meta_id	RefSeq	gene	$s_start	$s_end	.	+	.	ID=gene-S;Name=S;gbkey=Gene;gene=S;gene_biotype=protein_coding
$s_meta_id	RefSeq	CDS	$s_start	$s_end	.	+	0	ID=cds-S;Parent=gene-S;gbkey=CDS;product=nucleoprotein
EOF
            echo "  Created $S_META_GFF"
        else
            echo "Error: Could not extract CDS coordinates from $S_ORIGINAL_GFF"
        fi
    fi
else
    echo "Missing files needed for S segment metaconsensus GFF creation"
fi

# M segment: Extract CDS information and create metaconsensus GFF
if [ -f "$M_META" ] && [ -f "$M_ORIGINAL_GFF" ]; then
    # Get the sequence ID from metaconsensus FASTA
    m_meta_id=$(grep ">" "$M_META" | head -1 | sed 's/>//' | cut -d ' ' -f1)
    
    if [ -z "$m_meta_id" ]; then
        echo "Error: Could not extract sequence ID from $M_META"
    else
        echo "Creating GFF for M segment metaconsensus (ID: $m_meta_id)..."
        
        # Extract the start and end positions from the original GFF
        m_start=$(grep "CDS" "$M_ORIGINAL_GFF" | awk '{print $4}')
        m_end=$(grep "CDS" "$M_ORIGINAL_GFF" | awk '{print $5}')
        
        if [ -n "$m_start" ] && [ -n "$m_end" ]; then
            # Create a GFF file for the metaconsensus
            cat > "$M_META_GFF" << EOF
##gff-version 3
$m_meta_id	RefSeq	gene	$m_start	$m_end	.	+	.	ID=gene-M;Name=M;gbkey=Gene;gene=M;gene_biotype=protein_coding
$m_meta_id	RefSeq	CDS	$m_start	$m_end	.	+	0	ID=cds-M;Parent=gene-M;gbkey=CDS;product=glycoprotein_precursor
EOF
            echo "  Created $M_META_GFF"
        else
            echo "Error: Could not extract CDS coordinates from $M_ORIGINAL_GFF"
        fi
    fi
else
    echo "Missing files needed for M segment metaconsensus GFF creation"
fi

echo "Metaconsensus GFF files created." 