#!/usr/bin/env python3

## mito_close_consensus.py version 0.3 (development)
# author: Miles Benton
# created: 2024/02/12 09:18:03
# last modified: 2024/02/12 15:33:05

# this script takes input fasta from mito_closer.py and finds a potential consensus sequence to complete
# the mitochondrial genome and close the loop. Consensus and reverse consensus are written out to fasta,
# and sequence that differs significantly is printed to stdout (with parent readID) for inspection.

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import argparse

def get_consensus_sequence(extracted_sequences):
    forward_sequences = [seq['sequence'] for seq in extracted_sequences if seq['info'][2] == 'forward']
    reverse_sequences = [seq['sequence'] for seq in extracted_sequences if seq['info'][2] == 'reverse']
    
    # Calculate the length of each sequence
    forward_lengths = [len(seq) for seq in forward_sequences]
    reverse_lengths = [len(seq) for seq in reverse_sequences]

    # Find the most common length for forward and reverse sequences
    common_forward_length = Counter(forward_lengths).most_common(1)[0][0]
    common_reverse_length = Counter(reverse_lengths).most_common(1)[0][0]
    
    # Filter out sequences with lengths significantly different from the common lengths
    filtered_forward_sequences = [seq for seq in forward_sequences if len(seq) == common_forward_length]
    filtered_reverse_sequences = [seq for seq in reverse_sequences if len(seq) == common_reverse_length]

    # Retain sequences with lengths different from the common lengths
    filtered_out_sequences = [seq for seq in extracted_sequences if len(seq['sequence']) != common_forward_length and len(seq['sequence']) != common_reverse_length]

    # Reverse complement sequences in reverse orientation
    filtered_reverse_sequences = [str(Seq(seq).reverse_complement()) for seq in filtered_reverse_sequences]

    # Concatenate filtered sequences
    filtered_sequences = filtered_forward_sequences + filtered_reverse_sequences

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
    parser.add_argument("output_file", help="Path to the output FASTA file for consensus and reverse consensus sequences")
    args = parser.parse_args()

    # Read extracted sequences from the input FASTA file
    extracted_sequences = []
    with open(args.input_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            sequence_info = {'sequence': str(record.seq)}
            header_info = record.description.split(':')
            sequence_info['info'] = header_info
            extracted_sequences.append(sequence_info)

    # Generate consensus sequence and retain filtered out sequences
    consensus_sequence, filtered_out_sequences = get_consensus_sequence(extracted_sequences)
    
    # Compute reverse complement of consensus sequence
    consensus_sequence_rc = str(Seq(consensus_sequence).reverse_complement())

    # Write consensus and reverse complement sequences to a FASTA file
    consensus_record = SeqRecord(Seq(consensus_sequence), id="consensus", description="")
    reverse_consensus_record = SeqRecord(Seq(consensus_sequence_rc), id="reverse_consensus", description="")
    SeqIO.write([consensus_record, reverse_consensus_record], args.output_file, "fasta")
    
    print("Consensus and reverse complement sequences written to:", args.output_file)

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
