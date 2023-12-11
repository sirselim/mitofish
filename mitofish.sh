#!/bin/bash -l

## mitofish.sh version 0.3 (development)
# author: Miles Benton
# created: 2023/12/11 08:17:08
# last modified: 2023/12/12 09:56:56

# This script generates "baited" sequences from across a provided genome to search for them in fastq file(s)
# results are passed to stdout and can be redirected to fastq/fq or fastq.gz (via bgzip)

#TODO list
# - explore using flexiplex rather than seqkit? 

# depends on:
# seqkit

# Check if seqkit is installed
if ! command -v seqkit &> /dev/null; then
    echo "Error: 'seqkit' is not installed. Please install it before running this script."
    exit 1
fi

# version
VERSION="0.3dev"

# required input
threads="4"
genome_size=""
bait_length="60"
reference=""
fastq_input=""
output=""
mismatch="3"

# Function to display usage information
usage() {
    echo -e "Usage: $0 [-t, --threads <threads> (optional)] [-g, --genome-size <genome_size>] [-b, --bait-length <bait_length> (optional)] 
      [-m, --mismatch <mismatch>] [-r, --reference <reference>] [-f, --fastq-input <fastq_input>] [-o, --output <output>]\n"
    echo -e "$0 -h | --help will display the help menu\n"
    show_help
    exit 1
}

# Help function
show_help() {
    echo "mitofish - version: ${VERSION}"
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -h, --help          Show this help message and exit"
    echo "  -t, --threads       [optional] Number of CPU threads to use for processing, default: 4"
    echo "  -g, --genome-size   Approx size of the reference geneome being used, example: 16500 (for mammalian mt genomes)"
    echo "  -b, --bait-length   [optional] Length of bait sequence used for searching, default: 60"
    echo "  -m, --mismatch      [optional] Number of mismatches to allow in the bait sequences, default: 3"
    echo "  -r, --reference     Reference file to generate bait sequences"
    echo "  -f, --fastq-input   Path to fq|fastq file to extract reads from"
    echo "  -o, --output        File path for extracted reads to be written to"
}

# Process command-line arguments
while [ "$#" -gt 0 ]; do
    case "$1" in
        -h|--help)
            show_help
            exit 0
            ;;
        -t|--threads)
            shift
            threads="$1"
            ;;
        -g|--genome-size)
            shift
            genome_size="$1"
            ;;
        -b|--bait-length)
            shift
            bait_length="$1"
            ;;
        -m|--mismatch)
            shift
            mismatch="$1"
            ;;
        -r|--reference)
            shift
            reference="$1"
            ;;
        -f|--fastq-input)
            shift
            fastq_input="$1"
            ;;
        -o|--output)
            shift
            output="$1"
            ;;
        *)
            echo "Unknown option: $1"
            usage
            show_help
            exit 1
            ;;
    esac
    shift
done

# Check if required inputs are provided
if [ -z "$genome_size" ] || [ -z "$reference" ] || [ -z "$fastq_input" ] || [ -z "$output" ]; then
    echo -e "Error: All inputs are required."
    usage
fi

# displaying the user input
echo -e "\n... running mitofish ...\n"
echo -e "Selected parameters:"
echo -e "  - threads: $threads"
echo -e "  - genome-size: $genome_size"
echo -e "  - bait-length: $bait_length"
echo -e "  - mismatch: $mismatch"
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
seqkit grep -j "${threads}" -i -m "${mismatch}" \
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