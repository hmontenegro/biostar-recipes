set -ueo pipefail

# Obtain the parameters.
# Input reads.
INPUT={{reads.toc}}

# Estimated average fragment length.
FRAG_LEN={{fragment_length.value}}

# Estimated standard devitaion of fragment length.
FRAG_SD={{fragment_sd.value}}

# Reference transcriptome.
TRANSCRIPTS={{transcripts.value}}

# Directory with Kallisto results.
mkdir -p results

# Directory with intermediate files.
mkdir -p tmp

# Build the Kallisto index.
INDEX_DIR={{runtime.local_root}}/temp
mkdir -p ${INDEX_DIR}
INDEX=${INDEX_DIR}/{{transcripts.uid}}.idx


# Build the Kallisto index it does not already exist.
if [ ! -f ${INDEX} ]; then
    echo "Building the kallisto index."
    kallisto index -i ${INDEX} ${TRANSCRIPTS}
else
    echo "Found an existing kallisto index."
fi

# Calculate abundances using kallisto.
cat ${INPUT}| egrep "fastq|fq" | sort | parallel -j 4 kallisto quant -o results/{1/.}.out  -i ${INDEX} --single -l ${FRAG_LEN} -s ${FRAG_SD} {}

# Create a combined count table for all samples.

# Directory name is the sample name.
# Rename est_count column to sample name and extract it into a new file.
ls -d results/* | cut -f 2 -d "/" | parallel "sed -e 1s/est_counts/{.}/ results/{}/abundance.tsv | cut -f 4 >tmp/{}_counts.txt"

# Combine all counts into a single file.
paste tmp/*counts.txt > tmp/all.txt

# Add  transcript ids to the counts.
ls -d results/* | head -1 | parallel cat {}/abundance.tsv | cut -f 1 | paste - tmp/all.txt > sample_counts.txt
