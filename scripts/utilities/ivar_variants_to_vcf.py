#!/usr/bin/env python3

import sys
import os
import datetime

def read_fasta(fasta_file):
    """Read a FASTA file and return a dictionary with sequence names as keys and sequences as values."""
    sequences = {}
    current_name = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]  # Extract just the first word after '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_name:  # Add the last sequence
            sequences[current_name] = ''.join(current_seq)
    
    return sequences

def get_segment_id(segment_path):
    """Extract segment information from the path for SnpEff compatibility."""
    if "S_segment" in segment_path:
        return "JN232078.1"
    elif "M_segment" in segment_path:
        return "OR184986.1"
    elif "L_segment" in segment_path:
        return "JN232080.1"
    else:
        return None

def ivar_to_vcf(tsv_file, fasta_file):
    """Convert iVar TSV output to VCF format."""
    
    # Read reference sequence
    sequences = read_fasta(fasta_file)
    
    # Get reference ID 
    ref_id = list(sequences.keys())[0]  # Use first sequence as reference
    ref_len = len(sequences[ref_id])
    
    # Determine segment for SnpEff compatibility
    snpeff_id = get_segment_id(fasta_file)
    if snpeff_id:
        # If we found a segment ID, use it for SnpEff compatibility
        contig_id = snpeff_id
    else:
        # Otherwise use the original ID from FASTA
        contig_id = ref_id
    
    # Get today's date in YYYYMMDD format
    today = datetime.datetime.now().strftime("%Y%m%d")
    
    # Write VCF header
    print("##fileformat=VCFv4.2")
    print(f"##fileDate={today}")
    print("##source=iVar-to-VCF conversion script")
    print("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">")
    print("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">")
    print("##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand Bias\">")
    print("##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Counts for ref-forward, ref-reverse, alt-forward, and alt-reverse\">")
    print("##INFO=<ID=AA,Number=1,Type=String,Description=\"Amino acid change annotation\">")
    print("##INFO=<ID=CODON,Number=1,Type=String,Description=\"Codon change\">")
    print("##INFO=<ID=NONSYN,Number=0,Type=Flag,Description=\"Non-synonymous mutation\">")
    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    print("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">")
    print("##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">")
    print("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele Depth\">")
    print(f"##contig=<ID={contig_id},length={ref_len}>")
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE")
    
    # Process iVar TSV file
    with open(tsv_file, 'r') as f:
        # Skip header line
        next(f)
        
        for line in f:
            fields = line.strip().split('\t')
            
            # Skip if fewer than expected columns
            if len(fields) < 11:
                continue
            
            # Parsing values
            region = fields[0]
            pos = fields[1]
            ref = fields[2]
            alt = fields[3]
            ref_dp = fields[4]
            ref_rv = fields[5]
            ref_qual = fields[6]
            alt_dp = fields[7]
            alt_rv = fields[8]
            alt_qual = fields[9]
            alt_freq = fields[10]
            total_dp = fields[11] if len(fields) > 11 else str(int(ref_dp) + int(alt_dp))
            pass_filter = fields[13] if len(fields) > 13 else "TRUE"
            filter_status = "PASS" if pass_filter.upper() in ["TRUE", "PASS", "T", "YES", "Y"] else "FAIL"
            
            # Extract amino acid information if available
            aa_info = ""
            nonsyn_flag = ""
            if len(fields) > 14:
                feature = fields[14]
                if len(fields) > 18 and fields[16] != "NA" and fields[18] != "NA":
                    ref_aa = fields[16]
                    alt_aa = fields[18]
                    codon_change = f"{fields[15]}>{fields[17]}"
                    aa_info = f";AA={ref_aa}{fields[19]}{alt_aa};CODON={codon_change}"
                    if ref_aa != alt_aa:
                        nonsyn_flag = ";NONSYN"
            
            # Build INFO field
            info = f"DP={ref_dp};AF={alt_freq};DP4={ref_qual},{ref_dp},{alt_rv},{alt_qual}{aa_info}{nonsyn_flag}"
            
            # Build FORMAT and genotype fields
            format_field = "GT:DP:AF:AD"
            genotype = "0/1"  # Heterzygous genotype for all variants
            sample = f"{genotype}:{ref_dp}:{alt_freq}:{ref_dp},{alt_dp}"
            
            # Write VCF record using the SnpEff-compatible contig ID
            print(f"{contig_id}\t{pos}\t.\t{ref}\t{alt}\t.\t{filter_status}\t{info}\t{format_field}\t{sample}")

def main():
    # Check arguments
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <ivar_tsv> <reference_fasta>", file=sys.stderr)
        sys.exit(1)
    
    tsv_file = sys.argv[1]
    fasta_file = sys.argv[2]
    
    # Validate input files
    if not os.path.exists(tsv_file):
        print(f"Error: TSV file not found: {tsv_file}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(fasta_file):
        print(f"Error: FASTA file not found: {fasta_file}", file=sys.stderr)
        sys.exit(1)
    
    # Convert TSV to VCF
    ivar_to_vcf(tsv_file, fasta_file)

if __name__ == "__main__":
    main() 