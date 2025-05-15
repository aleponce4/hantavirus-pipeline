#!/bin/bash

# Create a simple GFF3 file suitable for iVar from the existing GFF information

# For S segment
S_SEGMENT_GFF="data/references/S_segment/JN232078.gff3"
S_SEGMENT_SIMPLE="data/references/S_segment/JN232078.simple.gff"

# For M segment
M_SEGMENT_GFF="data/references/M_segment/OR184986.gff3"
M_SEGMENT_SIMPLE="data/references/M_segment/OR184986.simple.gff"

# S segment: Extract CDS information and create simple GFF
if [ -f "$S_SEGMENT_GFF" ]; then
    # Extract CDS start and end positions (typically 43-1329 for S segment)
    s_start=$(grep -A1 "JN232078.1" "$S_SEGMENT_GFF" | grep "CDS" | awk '{print $4}')
    s_end=$(grep -A1 "JN232078.1" "$S_SEGMENT_GFF" | grep "CDS" | awk '{print $5}')
    s_product=$(grep -A1 "JN232078.1" "$S_SEGMENT_GFF" | grep "CDS" | grep -o "product=[^;]*" | cut -d= -f2)
    
    if [ -n "$s_start" ] && [ -n "$s_end" ]; then
        echo "Creating simple GFF for S segment: CDS from $s_start to $s_end"
        
        # Create a simple GFF file with the required format
        cat > "$S_SEGMENT_SIMPLE" << EOF
##gff-version 3
JN232078.1	RefSeq	gene	$s_start	$s_end	.	+	.	ID=gene-S;Name=S;gbkey=Gene;gene=S;gene_biotype=protein_coding
JN232078.1	RefSeq	CDS	$s_start	$s_end	.	+	0	ID=cds-S;Parent=gene-S;gbkey=CDS;product=$s_product
EOF
        echo "  Created $S_SEGMENT_SIMPLE"
    else
        echo "Could not extract CDS information for S segment"
    fi
else
    echo "S segment GFF file not found: $S_SEGMENT_GFF"
fi

# M segment: Extract CDS information and create simple GFF
if [ -f "$M_SEGMENT_GFF" ]; then
    # Extract CDS start and end positions (typically 52-3468 for M segment)
    m_start=$(grep -A1 "OR184986.1" "$M_SEGMENT_GFF" | grep "CDS" | awk '{print $4}')
    m_end=$(grep -A1 "OR184986.1" "$M_SEGMENT_GFF" | grep "CDS" | awk '{print $5}')
    m_product=$(grep -A1 "OR184986.1" "$M_SEGMENT_GFF" | grep "CDS" | grep -o "product=[^;]*" | cut -d= -f2)
    
    if [ -n "$m_start" ] && [ -n "$m_end" ]; then
        echo "Creating simple GFF for M segment: CDS from $m_start to $m_end"
        
        # Create a simple GFF file with the required format
        cat > "$M_SEGMENT_SIMPLE" << EOF
##gff-version 3
OR184986.1	RefSeq	gene	$m_start	$m_end	.	+	.	ID=gene-M;Name=M;gbkey=Gene;gene=M;gene_biotype=protein_coding
OR184986.1	RefSeq	CDS	$m_start	$m_end	.	+	0	ID=cds-M;Parent=gene-M;gbkey=CDS;product=$m_product
EOF
        echo "  Created $M_SEGMENT_SIMPLE"
    else
        echo "Could not extract CDS information for M segment"
    fi
else
    echo "M segment GFF file not found: $M_SEGMENT_GFF"
fi

echo "Simple GFF files created." 