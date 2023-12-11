#!/bin/bash

## mitofish.sh version 0.2
## author: Miles Benton
## created: 2023/12/11 08:17:08

# This script generates "baited" sequences from across a provided genome to search for them in fastq file(s)
# results are passed to stdout and can be redirected to fastq/fq or fastq.gz (via bgzip)

# depends on:
# seqkit

# Check if seqkit is installed
if ! command -v seqkit &> /dev/null; then
    echo "Error: 'seqkit' is not installed. Please install it before running this script."
    exit 1
fi

# required input
threads="4"
genome_size=""
bait_length="60"
reference=""
fastq_input=""
output=""

# Function to display usage information
usage() {
    echo "Usage: $0 [--threads <threads> (optional)] [--genome-size <genome_size>] [--bait-length <bait_length> (optional)] [--reference <reference>] [--fastq-input <fastq_input>] [--output <output>]"
    exit 1
}

# Process command-line arguments
while [ "$#" -gt 0 ]; do
    case "$1" in
        --threads)
            shift
            threads="$1"
            ;;
        --genome-size)
            shift
            genome_size="$1"
            ;;
        --bait-length)
            shift
            bait_length="$1"
            ;;
        --reference)
            shift
            reference="$1"
            ;;
        --fastq-input)
            shift
            fastq_input="$1"
            ;;
        --output)
            shift
            output="$1"
            ;;
        *)
            usage
            ;;
    esac
    shift
done

# Check if required inputs are provided
if [ -z "$genome_size" ] || [ -z "$reference" ] || [ -z "$fastq_input" ] || [ -z "$output" ]; then
    echo "Error: All inputs are required."
    usage
fi

# displaying the user input
echo -e "\n... running mitofish ...\n"
echo -e "Selected parameters:"
echo -e "  - threads: $threads"
echo -e "  - genome-size: $genome_size"
echo -e "  - bait-length: $bait_length"
echo -e "  - reference: $reference"
echo -e "  - fastq-input: $fastq_input"
echo -e "  - output: $output"

# defined values
START=1  # where to start the first interval
END="${genome_size}"  # adjust the 'end' value as needed, this is roughly mt genome size
# define spacing of baits across reference
SPACING=$(($genome_size / 16))
echo -e "  - spacing of bait sequences has been calculated to: $SPACING bp\n"   # how spaced apart the seqences should be
INTERVAL_LENGTH="${bait_length}"  # this is the length of the sequence we want to grep for


# create an empty list to populate with intervals
INTERVAL_LIST=()

# loop to generate intervals
for ((INTERVAL_START = START; INTERVAL_START <= END; INTERVAL_START += SPACING)); do
  INTERVAL_END=$((INTERVAL_START + INTERVAL_LENGTH - 1)) 
  # ensure that the interval doesn't exceed the 'end' value
  if ((INTERVAL_END > END)); then
    INTERVAL_END=$END
  fi
  INTERVAL="$INTERVAL_START-$INTERVAL_END"
  INTERVAL_LIST+=("${INTERVAL}")
done

# echo "(${INTERVAL_LIST[*]})"  # debugging

## create a list of sequences to bait out mt reads from fastq's
# create an empty list to populate with sequences
SEQ_LIST=()

# generate the sub-seqs
for INT in "${INTERVAL_LIST[@]}"; do
  CURRENT_SEQ="$(tail -n +2 "${reference}" | cut -c "${INT}")";
  SEQ_LIST+=("${CURRENT_SEQ}");
done

# echo "(${SEQ_LIST[*]})"  # debugging

## there is a issue where the last position of the array can contain a sequence that
## is not of the defined length. This below section should check for this and remove
## the offending sequence if it's the case.

# Threshold value
threshold="${bait_length}"

# Create a new array to store filtered values
filtered_array=()

# Loop through the elements of the original array
for value in "${SEQ_LIST[@]}"; do
    # Check if the value is greater than or equal to the threshold
    if [ "${#value}" -ge "${threshold}" ]; then
        # Add the value to the filtered array
        filtered_array+=("${value}")
    fi
done

### DEBUGGING
# # Display the original and filtered arrays
# echo -e "\n"
# echo -e "Original Array: ${SEQ_LIST[@]}"
# echo -e "Filtered Array (values greater than or equal to $threshold): ${filtered_array[@]}"

# replace original array with filtered
SEQ_LIST="${filtered_array}"

## use seqkit grep to pull out mt reads based on a prior reference or near-related species
echo -e "... seqkit running ..."
# seqkit grep for mt reads
seqkit grep -j "${threads}" -i -m 3 \
  -p "${SEQ_LIST[0]}" \
  -p "${SEQ_LIST[1]}" \
  -p "${SEQ_LIST[2]}" \
  -p "${SEQ_LIST[3]}" \
  -p "${SEQ_LIST[4]}" \
  -p "${SEQ_LIST[5]}" \
  -p "${SEQ_LIST[6]}" \
  -p "${SEQ_LIST[7]}" \
  -p "${SEQ_LIST[8]}" \
  -p "${SEQ_LIST[9]}" \
  -p "${SEQ_LIST[10]}" \
  -p "${SEQ_LIST[11]}" \
  -p "${SEQ_LIST[12]}" \
  -p "${SEQ_LIST[13]}" \
  -p "${SEQ_LIST[14]}" \
  -p "${SEQ_LIST[15]}" \
  "${fastq_input}" --out-file "${output}"

# finished
# report some stats
echo -e "\n"
echo -e "... seqkit statistics ..."
seqkit stats "${output}"
READS=$(seqkit stats ${output} | awk '{print $4}' | tail -n 1)
echo -e "\nmitofish output ${READS} potential mitochondrial reads to: ${output}"

# close
echo -e "\n... mitofish completed ...\n"
##/END