#!/usr/bin/env python3

## mito_close_consensus.py version 0.1 (development)
# author: Miles Benton
# created: 2024/02/12 09:18:03
# last modified: 2024/02/12 12:03:27

from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
import argparse

def get_consensus_sequence(extracted_sequences):
    # Calculate the length of each sequence
    sequence_lengths = [len(seq) for seq in extracted_sequences]
    
    # Find the most common length
    common_length = Counter(sequence_lengths).most_common(1)[0][0]
    
    # Filter out sequences with lengths significantly different from the common length
    filtered_sequences = [seq for seq in extracted_sequences if len(seq) == common_length]
    filtered_out_sequences = [seq for seq in extracted_sequences if len(seq) != common_length]
    
    # Determine the predominant orientation (forward or reverse)
    orientations = [seq.isupper() for seq in filtered_sequences]
    forward_count = orientations.count(True)
    reverse_count = orientations.count(False)
    predominant_orientation = "forward" if forward_count > reverse_count else "reverse"

    # Reverse complement sequences in reverse orientation
    if predominant_orientation == "reverse":
        filtered_sequences = [str(Seq(seq).reverse_complement()) for seq in filtered_sequences]

    # Compare sequences and determine consensus
    consensus_sequence = ""
    for positions in zip(*filtered_sequences):
        counts = Counter(positions)
        consensus_base = counts.most_common(1)[0][0]
        consensus_sequence += consensus_base

    return consensus_sequence, filtered_out_sequences

def main():
    parser = argparse.ArgumentParser(description="Generate consensus sequence from extracted sequences in a FASTA file")
    parser.add_argument("input_file", help="Path to the FASTA file containing extracted sequences")
    args = parser.parse_args()

    # Read extracted sequences from the input FASTA file
    extracted_sequences = []
    with open(args.input_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            extracted_sequences.append(str(record.seq))

    # Generate consensus sequence
    consensus_sequence, filtered_out_sequences = get_consensus_sequence(extracted_sequences)
    
    # Print consensus sequence
    print("Consensus sequence:", consensus_sequence)

    # Compute and print reverse complement of consensus sequence
    consensus_sequence_rc = str(Seq(consensus_sequence).reverse_complement())
    print("Reverse complement of consensus sequence:", consensus_sequence_rc)
    
    # Print filtered out sequences
    print("\nFiltered out sequences (due to different lengths):")
    for seq in filtered_out_sequences:
        print(seq)

if __name__ == "__main__":
    main()



