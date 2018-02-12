set -ueo pipefail

# Input reads.
INPUT={{reads.toc}}

# Reference Transcriptome.
TRANSCRIPTS={{transcripts.value}}

# Library type.
LIBRARY={{library.value}}

# Tool to be used.
TOOL={{tool.value}}

# Fraction to sub-sample the data.
FRACTION={{fraction.value}}

# Estimated average fragment length.
FRAG_LEN={{fragment_length.value}}

# Estimated standard deviation of fragment length.
FRAG_SD={{fragment_sd.value}}

# Run log to redirect unwanted output.
RUNLOG=runlog/runlog.txt

# This directory should already exist.
mkdir -p runlog

# Directory with results.
mkdir -p results

# Directory with intermediate files.
mkdir -p tmp

# Sub-sample data.
if [ ${FRACTION} != 1 ]; then

    # Directory with sub-sampled data
    READS=reads
    mkdir -p $READS

    # Subsample data.
    echo "Randomly sampling $FRACTION fraction of data"
    cat ${INPUT} | egrep "fastq|fq" | sort |  parallel "seqtk sample -s 11 {} $FRACTION >$READS/{/.}.fq"

    # The name of the new table of contents.
    INPUT=$READS/toc.txt

    # Create table of contents with sub-sampled data.
    ls -1 $READS/*.fq > $INPUT
fi

# Index directory.
INDEX_DIR={{runtime.local_root}}/indices
mkdir -p ${INDEX_DIR}

# Kallisto index.
KINDEX=${INDEX_DIR}/{{transcripts.uid}}_kallisto.idx

# Salmon index.
SINDEX=${INDEX_DIR}/{{transcripts.uid}}_salmon.idx

# Build the Kallisto index if needed.
if [ "$TOOL" == "kallisto" ] && [ ! -f "$KINDEX" ]; then
    echo "Building the kallisto index."
    kallisto index -i ${SINDEX} ${TRANSCRIPTS} >> $RUNLOG 2>&1
fi

# Build the Salmon index if needed.
if [ "$TOOL" == "salmon" ] && [ ! -f "$SINDEX" ]; then
    echo "Building the salmon index."
    salmon index -i ${SINDEX} -t ${TRANSCRIPTS} -p 2 >> $RUNLOG 2>&1
fi

# Kallisto in single end mode.
if [ "$TOOL" == "kallisto" ] && [ "$LIBRARY" == "SE" ]; then
    echo  "Estimating transcript abundances using $TOOL in single end mode."
    cat ${INPUT}| egrep "fastq|fq" | sort | parallel -j 4 kallisto quant -o results/{1/.}.out  -i ${KINDEX} --single
    -l ${FRAG_LEN} -s ${FRAG_SD} {} >> $RUNLOG 2>&1
fi

# Kallisto in paired end mode.
if [ "$TOOL" == "kallisto" ] && [ "$LIBRARY" == "PE" ]; then
    echo "Estimating transcript abundances usinf $TOOL in paired end mode."
    cat ${INPUT}| egrep "fastq|fq" | sort | parallel -N 2  -j 4 kallisto quant -o results/{1/.}.out  -i ${KINDEX} {1}
     {2} >> $RUNLOG 2>&1
fi

# Salmon in single end mode.
if [ "$TOOL" == "salmon" ] && [ "$LIBRARY" == "SE" ]; then
    echo  "Estimating transcript abundances using $TOOL in single end mode."
    cat ${INPUT}| egrep "fastq|fq" | sort | parallel -j 4 salmon quant -o results/{1/.}.out  -i ${SINDEX} --libType A
     -r {} >> $RUNLOG 2>&1
fi


# Salmon in paired end mode.
if [ "$TOOL" == "salmon" ] && [ "$LIBRARY" == "PE" ]; then
    echo  "Estimating transcript abundances using $TOOL in single end mode."
    cat ${INPUT}| egrep "fastq|fq" | sort | parallel -N 2 -j 4 salmon quant -o results/{1/.}.out  -i ${SINDEX}
    --libType A -1 {1} -2 {2} >> $RUNLOG 2>&1
fi


# Create a combined count table for all samples.

if [ "$TOOL" == "kallisto" ]; then
    # Directory name is the sample name.
    # Rename est_count column to sample name and extract it into a new file.
    ls -d results/* | cut -f 2 -d "/" | parallel "sed -e 1s/est_counts/{.}/ results/{}/abundance.tsv | cut -f 4 >tmp/{}_counts.txt"
fi

if [ "$TOOL" == "salmon" ]; then
    # Directory name is the sample name.
    # Rename est_count column to sample name and extract it into a new file.
    ls -d results/* | cut -f 2 -d "/" | parallel "sed -e 1s/NumReads/{.}/ results/{}/quant.sf | cut -f 5 >tmp/{}_counts.txt"
fi

# Combine all counts into a single file.
paste tmp/*counts.txt > tmp/all.txt

# Add  transcript ids to the counts.
ls -d results/* | head -1 | parallel cat {}/abundance.tsv | cut -f 1 | paste - tmp/all.txt > combined_abundance.txt



