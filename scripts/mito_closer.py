#!/usr/bin/env python3

## mito_closer.py version 0.1 (development)
# author: Miles Benton
# created: 2024/02/12 09:18:03
# last modified: 2024/02/12 12:02:34

# This script takes a mito assembly as input (fasta) and then uses 30bp from either end
# to try and locate sequence to help complete the genome and close the circle. It relies
# on the logic that if the mt is circular then if you take a selection from either end
# of a linear assembly they will either touch (nothing between them), or there will be
# an amount of sequence between them. This allows the true length to be determined and 
# essentially closes the genome.

#TODO list
# - provide proper help options
# - add an option to toggle verbose
# - test more!!

from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import gzip

def extract_target_sequences(fasta_file, fastq_file, output_file):
    # Read the 30bp sequences from the start and end of the FASTA file
    with open(fasta_file, "r") as f:
        fasta_sequences = SeqIO.parse(f, "fasta")
        for record in fasta_sequences:
            start_sequence = str(record.seq[:30])
            end_sequence = str(record.seq[-30:])
            break  # Assuming only one sequence in the FASTA file

    # Find reverse complements of start and end sequences
    start_sequence_rc = str(Seq(start_sequence).reverse_complement())
    end_sequence_rc = str(Seq(end_sequence).reverse_complement())

    print("... running mito closer ...\n")
    print("Start sequence:", start_sequence)
    print("End sequence:", end_sequence)
    print("Start sequence RC:", start_sequence_rc)
    print("End sequence RC:", end_sequence_rc)
    print("\n")

    # Create a dictionary to store extracted sequences with their corresponding read IDs
    extracted_sequences = {}

    # Search for reads containing either forward or reverse complement of start and end sequences in the FASTQ file
    with gzip.open(fastq_file, "rt") as f:
        fastq_sequences = SeqIO.parse(f, "fastq")
        for record in fastq_sequences:
            sequence = str(record.seq)
            forward_match = (start_sequence in sequence and end_sequence in sequence)
            reverse_match = (start_sequence_rc in sequence and end_sequence_rc in sequence)
            if forward_match or reverse_match:
                print("Read ID:", record.id)
                if forward_match:
                    start_index = sequence.index(start_sequence)
                    end_index = sequence.index(end_sequence)
                else:  # reverse_match
                    start_index = sequence.index(start_sequence_rc)
                    end_index = sequence.index(end_sequence_rc)
                
                print("Start index:", start_index)
                print("End index:", end_index)
                
                if start_index > end_index:
                    start_index, end_index = end_index, start_index
                
                extracted_seq = sequence[start_index+len(start_sequence):end_index]
                print("Extracted sequence:", extracted_seq)
                
                if extracted_seq:
                    extracted_sequences[record.id] = extracted_seq

    # Write the extracted sequences to the output file in FASTA format
    with open(output_file, "w") as f:
        for read_id, seq in extracted_sequences.items():
            f.write(f">parent_read:{read_id}\n{seq}\n")
            
    # Print information about the number of extracted sequences and the output file path
    print("\nNumber of extracted sequences:", len(extracted_sequences))
    print("Extracted sequences written to:", output_file)

def main():
    parser = argparse.ArgumentParser(description="Extract sequences between start and end sequences from FASTQ file")
    parser.add_argument("fasta_file", help="Path to the FASTA file containing start and end sequences")
    parser.add_argument("fastq_file", help="Path to the FASTQ file (gzip-compressed) to search for sequences")
    parser.add_argument("output_file", help="Path to the output file to store extracted sequences")
    args = parser.parse_args()

    extract_target_sequences(args.fasta_file, args.fastq_file, args.output_file)

if __name__ == "__main__":
    main()
