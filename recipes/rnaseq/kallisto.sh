set -ueo pipefail

# Input reads.
INPUT={{reads.toc}}

# Reference Transcriptome.
TRANSCRIPTS={{transcripts.value}}

# Library type.
LIBRARY={{library.value}}

# Fraction to sub-sample the data.
FRACTION={{sampling.value}}

# Sub-sample data.
if [ ${FRACTION} == 1 ]; then
    # DATA and $INPUT are same if all reads are selected.
    DATA=$INPUT
else
    # Directory with sub-sampled data.
    SAMPLED=sampled
    mkdir -p $SAMPLED

    # Table of contents with sub-sampled data.
    DATA=$SAMPLED/sampled_toc.txt

    # Subsample data.
    echo "Randomly sampling $FRACTION fraction of data"
    cat ${INPUT} | egrep "fastq|fq" | sort |  parallel "seqtk sample -s 11 {} $FRACTION >$SAMPLED/{/.}_sampled$FRACTION.fq"

    # Create table of contents with sub-sampled data.
    ls `pwd`/$SAMPLED/*.fq >$DATA
fi

# Kallisto index directory.
INDEX_DIR={{runtime.local_root}}/indices
mkdir -p ${INDEX_DIR}

# Kallisto index.
INDEX=${INDEX_DIR}/{{transcripts.uid}}.idx

# Build the Kallisto index if it does not already exist.
if [ ! -f ${INDEX} ]; then
    echo "Building the kallisto index."
    kallisto index -i ${INDEX} ${TRANSCRIPTS}
else
    echo "Found an existing kallisto index."
fi

# Directory with Kallisto results.
mkdir -p results

# Directory with intermediate files.
mkdir -p tmp

# Calculate abundances using kallisto.
if [ ${LIBRARY} == "PE" ]; then
    echo "Running kallisto quant algorithm in paired end mode."
    cat ${DATA}| egrep "fastq|fq" | sort | parallel -N 2  -j 4 kallisto quant -o results/{1/.}.out  -i ${INDEX} {1} {2}
else
    # Obtain additional parameters for SE.
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



