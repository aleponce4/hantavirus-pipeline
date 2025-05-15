#!/usr/bin/env python3
"""
Debug script to check primer binding
"""
from Bio import SeqIO
from Bio.Seq import Seq
import os

# Load the M segment consensus sequence
consensus_path = "data/references/M_segment/metaconsensus.fasta"
consensus_seq = str(SeqIO.read(consensus_path, "fasta").seq)

# MF2864 primer details from CSV
primer_name = "MF2864"
primer_seq = "CTA CAA CCT TTG AGG CCC AGA CAA AG"
primer_region = "341-366F"
csv_overall_region = "2865-3302"

# Clean the sequence
primer_seq_clean = ''.join(c for c in primer_seq if c.upper() in 'ACGT')

# Check both positions
position_341 = consensus_seq[340:365]  # 0-based, checking where it says it should be
position_2864 = consensus_seq[2863:2863+len(primer_seq_clean)]  # 0-based, checking where name implies it should be

# Find best match by sliding window
best_match_pos = -1
best_match_score = 0

for i in range(len(consensus_seq) - len(primer_seq_clean) + 1):
    ref_segment = consensus_seq[i:i+len(primer_seq_clean)]
    match_count = sum(a == b for a, b in zip(primer_seq_clean, ref_segment))
    match_score = match_count / len(primer_seq_clean)
    
    if match_score > best_match_score:
        best_match_score = match_score
        best_match_pos = i

# Print results
print(f"Debugging primer: {primer_name} ({primer_region})")
print(f"Primer sequence: {primer_seq_clean}")
print(f"Overall region in CSV: {csv_overall_region}")
print(f"\nPosition 341-366 match check:")
print(f"Sequence at positions 341-366: {position_341}")
print(f"Match score: {sum(a == b for a, b in zip(primer_seq_clean, position_341)) / len(primer_seq_clean):.2f}")

print(f"\nPosition 2864 match check:")
print(f"Sequence at positions 2864-{2864+len(primer_seq_clean)-1}: {position_2864}")
print(f"Match score: {sum(a == b for a, b in zip(primer_seq_clean, position_2864)) / len(primer_seq_clean):.2f}")

print(f"\nBest match found:")
print(f"Best position: {best_match_pos+1}-{best_match_pos+len(primer_seq_clean)}")
print(f"Best match score: {best_match_score:.2f}")
print(f"Sequence at best match: {consensus_seq[best_match_pos:best_match_pos+len(primer_seq_clean)]}")

# Print alignment of primer with best match
print("\nAlignment:")
print(f"Primer: {primer_seq_clean}")
matches = ""
for a, b in zip(primer_seq_clean, consensus_seq[best_match_pos:best_match_pos+len(primer_seq_clean)]):
    matches += "|" if a == b else " "
print(f"Match:  {matches}")
print(f"Target: {consensus_seq[best_match_pos:best_match_pos+len(primer_seq_clean)]}")

# Check reverse primer
rev_name = "MR3274"
rev_seq = "GGCAGAATATGTGGAGTCAGTCAAATTAAATG"
rev_region = "3274-3306R"
rev_seq_clean = ''.join(c for c in rev_seq if c.upper() in 'ACGT')
rev_comp = str(Seq(rev_seq_clean).reverse_complement())

# Find best match for reverse primer
best_match_pos_rev = -1
best_match_score_rev = 0

for i in range(len(consensus_seq) - len(rev_comp) + 1):
    ref_segment = consensus_seq[i:i+len(rev_comp)]
    match_count = sum(a == b for a, b in zip(rev_comp, ref_segment))
    match_score = match_count / len(rev_comp)
    
    if match_score > best_match_score_rev:
        best_match_score_rev = match_score
        best_match_pos_rev = i

print(f"\n\nDebugging reverse primer: {rev_name} ({rev_region})")
print(f"Reverse primer sequence: {rev_seq_clean}")
print(f"Reverse complement: {rev_comp}")
print(f"Best match position: {best_match_pos_rev+1}-{best_match_pos_rev+len(rev_comp)}")
print(f"Best match score: {best_match_score_rev:.2f}")

# Calculate amplicon size based on best binding sites
if best_match_pos >= 0 and best_match_pos_rev >= 0:
    amplicon_size_from_binding = best_match_pos_rev + len(rev_comp) - best_match_pos
    print(f"\nAmplicon size if using best binding sites: {amplicon_size_from_binding} bp")

# Calculate amplicon size based on CSV positions
amplicon_size_from_csv = 3306 - 341 + 1
print(f"Amplicon size from CSV positions (341 to 3306): {amplicon_size_from_csv} bp")

# Calculate amplicon size based on primer names
amplicon_size_from_names = 3274 - 2864 + 33  # Adding reverse primer length
print(f"Amplicon size if positions matched names (2864 to 3306): {amplicon_size_from_names} bp") 