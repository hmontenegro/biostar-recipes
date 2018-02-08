set -ueo pipefail

# Input reads.
INPUT={{reads.toc}}

# Reference Transcriptome.
TRANSCRIPTS={{transcripts.value}}

# Protocol.
PROTOCOL={{protocol.value}}

# Fraction to sub-sample.
FRACTION={{sub_sample.value}}

# Directory with intermediate files.
mkdir -p tmp

# Sub-sample data.
if [ ${FRACTION} == 1.0 ]; then
    # DATA and $INPUT are same if $FRACTION=100%
    DATA=$INPUT
else
    # Table of contents with sub-sampled data.
    DATA=tmp/subset_toc.txt

    # Subsample data.
    cat ${INPUT} | egrep "fastq|fq" | sort |  parallel "seqtk sample -s 11 {} $FRACTION >tmp/{/.}_subset$FRACTION.fq"

    # Create table of contents with sub-sampled data.
    ls `pwd`/tmp/*.fq >$DATA
fi

# Kallisto index directory.
INDEX_DIR={{runtime.local_root}}/temp
mkdir -p ${INDEX_DIR}

# Kallisto index.
INDEX=${INDEX_DIR}/{{transcripts.uid}}.idx

# Build the Kallisto index it does not already exist.

if [ ! -f ${INDEX} ]; then
    echo "Building the kallisto index."
    kallisto index -i ${INDEX} ${TRANSCRIPTS}
else
    echo "Found an existing kallisto index."
fi

# Directory with Kallisto results.
mkdir -p results

# Calculate abundances using kallisto.
if [ ${PROTOCOL} == "paired" ]; then
    echo "Running kallisto quant algorithm in paired end mode."
    cat ${DATA}| egrep "fastq|fq" | sort | parallel -N 2  -j 4 kallisto quant -o results/{1/.}.out  -i ${INDEX} {1} {2}
else
    # Obtain additional parameters.
    # Estimated average fragment length.
    FRAG_LEN={{fragment_length.value}}

    # Estimated standard devitaion of fragment length.
    FRAG_SD={{fragment_sd.value}}

    # Run kallisto.
    echo "Running kallisto quant algorithm in single end mode."
    cat ${DATA}| egrep "fastq|fq" | sort | parallel -j 4 kallisto quant -o results/{1/.}.out  -i ${INDEX} --single -l ${FRAG_LEN} -s ${FRAG_SD} {}
fi

# Create a combined count table for all samples.

# Directory name is the sample name.
# Rename est_count column to sample name and extract it into a new file.
ls -d results/* | cut -f 2 -d "/" | parallel "sed -e 1s/est_counts/{.}/ results/{}/abundance.tsv | cut -f 4 >tmp/{}_counts.txt"

# Combine all counts into a single file.
paste tmp/*counts.txt > tmp/all.txt

# Add  transcript ids to the counts.
ls -d results/* | head -1 | parallel cat {}/abundance.tsv | cut -f 1 | paste - tmp/all.txt > sample_counts.txt



