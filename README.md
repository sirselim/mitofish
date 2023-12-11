# mitofish
 This script generates "baited" sequences from across a provided genome to search for them in fastq file(s)


## usage

```sh
./mitofish.sh --genome-size 49300 \
  --reference ./test_data/test_reference.fasta \
  --fastq-input ./test_data/test_reads.fastq.gz \
  --threads 8 \
  --bait-length 30 \
  --output mitofish_results.fastq.gz
```

```sh
$ ./mitofish.sh --help
mitofish - version: 0.3dev
Usage: ./mitofish.sh [options]
Options:
  -h, --help          Show this help message and exit
  -t, --threads       [optional] Number of CPU threads to use for processing, default: 4
  -g, --genome-size   Approx size of the reference geneome being used, example: 16500 (for mammalian mt genomes)
  -b, --bait-length   [optional] Length of bait sequence used for searching, default: 60
  -m, --mismatch      [optional] Number of mismatches to allow in the bait sequences, default: 3
  -r, --reference     Reference file to generate bait sequences
  -f, --fastq-input   Path to fq|fastq file to extract reads from
  -o, --output        File path for extracted reads to be written to
```