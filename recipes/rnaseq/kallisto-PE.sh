set -ueo pipefail

# Obtain the parameters.
INPUT={{reads.toc}}
TRANSCRIPTS={{transcripts.value}}

# Directory with Kallisto results.
mkdir -p results

# Build the Kallisto index.
INDEX_DIR={{runtime.local_root}}/temp
mkdir -p ${INDEX_DIR}
INDEX=${INDEX_DIR}/{{transcripts.uid}}


# Build the Kallisto index it does not already exist.
if [ ! -f "$INDEX.idx" ]; then
    echo "Building the kallisto index."
    kallisto index -i ${INDEX} ${TRANSCRIPTS}
else
    echo "Found an existing kallisto index."
fi

# Calculate abundances using kallisto.
cat ${INPUT}| egrep "fastq|fq" | sort | parallel -N 2 -j 4 kallisto quant -o results/{1/.}.out  -i ${INDEX} {1} {2}

# Create a combined abundance table for all samples.

# Get all counts with est_count column renamed with sample_name
ls -d results/* | cut -f 2 -d "/" |parallel "sed -e 1s/est_counts/{.}/ results/{}/abundance.tsv | cut -f 4 >{}_counts.txt"

# Combine counts from all samples.
paste *counts.txt >all.txt

# Add  transcripts ids to the counts.
ls -d results/* | head -1 | parallel cat {}/abundance.tsv | cut -f 1 | paste - all.txt > raw_abundance.txt

# Remove intermediate files.
rm -f all.txt *counts.txt


