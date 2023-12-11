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
