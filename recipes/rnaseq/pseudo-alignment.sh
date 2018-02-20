set -ueo pipefail

# Input reads.
INPUT={{reads.toc}}

# Reference Transcriptome.
TRANSCRIPTS={{transcripts.value}}

# Library protocol.
LIBRARY={{library.value}}

# Library type
LIBRARY_TYPE={{library_tpe.value}}

# PE read orientation
ORIENTATION={{orientation.value}}

# Tool to be used.
TOOL={{tool.value}}

# Fraction to sub-sample the data.
FRACTION={{fraction.value}}

# Estimated average fragment length.
FRAG_LEN={{fragment_length.value}}

# Estimated standard deviation of fragment length.
FRAG_SD={{fragment_sd.value}}

# Run log to redirect unwanted output.
RUNLOG=runlog/runlog.

# This directory should already exist.
mkdir -p runlog

# Directory with results.
mkdir -p results

# Directory with intermediate files.
mkdir -p tmp

# Directory with sub-sampled data.
READS=reads
mkdir -p $READS

# The name of the input files.
FILES=${READS}/files.txt

#Create file list
cat ${INPUT} | egrep 'fastq|fq'|sort >${FILES}

# Sub-sample data.
if [ ${FRACTION} != 1 ]; then

    # The random seed for sampling.
    SEED=$((1 + RANDOM % 1000))

    # Generate a random sample of each input file.
    echo "Sampling fraction=$FRACTION of data with random seed=$SEED"
    cat ${FILES} |  parallel "seqtk sample -s $SEED {} $FRACTION >$READS/{/.}.fq"

    # Create table of contents with sub-sampled data.
    ls -1 $READS/*.fq > ${FILES}
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

    if [ "$LIBRARY_TYPE" == "unstranded" ]; then
        cat ${INPUT}| egrep "fastq|fq" | sort | parallel -N 2  -j 4 kallisto quant -o results/{1/.}.out  -i ${KINDEX} {1}
     {2} >> $RUNLOG 2>&1

    else
        # run kallisto in stranded mode.
        cat ${INPUT}| egrep "fastq|fq" | sort | parallel -N 2  -j 4 kallisto quant -o results/{1/.}.out --${LIBRARY_TYPE} -i ${KINDEX} {1}
     {2} >> $RUNLOG 2>&1
     fi

fi

# Salmon in single end mode.
if [ "$TOOL" == "salmon" ] && [ "$LIBRARY" == "SE" ]; then

    echo  "Estimating transcript abundances using $TOOL in single end mode."
    cat ${INPUT}| egrep "fastq|fq" | sort | parallel -j 4 salmon quant -o results/{1/.}.out  -i ${SINDEX} --libType ${LIBRARY_TYPE}
     -r {} >> $RUNLOG 2>&1
fi


# Salmon in paired end mode.
if [ "$TOOL" == "salmon" ] && [ "$LIBRARY" == "PE" ]; then

    echo  "Estimating transcript abundances using $TOOL in single end mode."

    if ["${LIB_TYPE}" == "fr-stranded" ]; then
        LIB_TYPE_STR=${ORIENTATION}"SF"
    fi

    if ["${LIB_TYPE}" == "rf-stranded" ]; then
        LIB_TYPE_STR=${ORIENTATION}"SR"
    fi

    cat ${INPUT}| egrep "fastq|fq" | sort | parallel -N 2 -j 4 salmon quant -o results/{1/.}.out  -i ${SINDEX}
    --libType ${LIB_TYPE_STR} -1 {1} -2 {2} >> $RUNLOG 2>&1

fi


# Create a combined count table for all samples.

if [ "$TOOL" == "kallisto" ]; then
    # Directory name is the sample name.
    # Rename est_count column to sample name and extract it into a new file.
    ls -d results/* | cut -f 2 -d "/" | parallel "sed -e 1s/est_counts/{.}/ results/{}/abundance.tsv | cut -f 4 >tmp/{}_counts.txt"
fi

if [ "$TOOL" == "salmon" ]; then
    # Directory name is the sample name.
    # Rename NumReads column to sample name and extract it into a new file.
    ls -d results/* | cut -f 2 -d "/" | parallel "sed -e 1s/NumReads/{.}/ results/{}/quant.sf | cut -f 5 >tmp/{}_counts.txt"
fi

# Combine all counts into a single file.
paste tmp/*counts.txt > tmp/all.txt

# Add  transcript ids to the counts.
ls -d results/* | head -1 | parallel cat {}/abundance.tsv | cut -f 1 | paste - tmp/all.txt > combined_abundance.txt



