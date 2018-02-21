set -ueo pipefail

# Obtain the parameters.
INPUT={{reads.toc}}
GENOME={{genome.value}}
LIBRARY={{library.value}}
ALIGNER={{aligner.value}}
FRACTION={{fraction.value}}

# The first time around users must wait until this index is created
# before starting another job that uses this same index.
INDEX_DIR={{runtime.local_root}}/indices
mkdir -p ${INDEX_DIR}
INDEX=${INDEX_DIR}/{{genome.uid}}

# This directory will store the alignment.
mkdir -p bam

# Directory with sub-sampled data.
READS=reads
mkdir -p $READS

# The name of the input files.
FILES=${READS}/files.txt

# Run log to redirect unwanted output.
RUNLOG=runlog/runlog.txt

# Number of processes.
PROC=4

# This directory should already exist.
mkdir -p runlog

# Generate the input file names.
cat ${INPUT} | sort | egrep "fastq|fq" > ${FILES}

# This recipe uses the Django templating engine conditionals to render the script.
{% if fraction.value != "1" %}

    # Sub-sampling step.

    # The random seed for sampling.
    SEED=$((1 + RANDOM % 1000))

    # Generate a random sample of each input file.
    echo "Sampling fraction=$FRACTION of data with random seed=$SEED"
    cat ${FILES} | parallel -j $PROC "seqtk sample -s $SEED {} $FRACTION > $READS/{/.}.fq"

    # Create a new table of contents with sub-sampled data.
    ls -1 $READS/*.fq > ${FILES}

{% endif %}

{# Generate script depending on the  aligner #}

{% if aligner.value == "bwa" %}

    # Build the BWA index if needed.
    if [ ! -f "$INDEX.bwt" ]; then
        echo "Building the bwa index."
        bwa index -p ${INDEX} ${GENOME} >> $RUNLOG 2>&1
    fi

    {% if library.value == "SE" %}
        # Run bwa in single end mode.
        cat ${FILES} | parallel -j $PROC "bwa mem ${INDEX} {1} 2>> $RUNLOG | samtools sort > bam/{1/.}.bam"
    {% else %}
        # Run bwa in paired end mode.
        cat ${FILES} | parallel -N 2 -j $PROC "bwa mem ${INDEX} {1} {2} 2>> $RUNLOG | samtools sort > bam/{1/.}.bam"
    {% endif %}

{% endif %}

{% if aligner.value == "bowtie2" %}

    # Build the Bowtie2 index if needed.
    if [ ! -f "$INDEX.bt2" ]; then
        echo "Building the bowtie2 index."
        bowtie2-build  ${GENOME} ${INDEX} >> $RUNLOG 2>&1
    fi

    {% if library.value == "SE" %}
        # Run bowtie2 in single end mode.
        cat ${FILES} | parallel -j $PROC "bowtie2 -x ${INDEX} --sensitive-local -U {1} 2>> $RUNLOG | samtools sort > bam/{1/.}.bam"
    {% else %}
        # Run bowtie2 in paired end mode.
        cat ${FILES} | parallel -N 2 -j $PROC "bowtie2 -x ${INDEX} --sensitive-local -1 {1} -2 {2} 2>> $RUNLOG | samtools sort > bam/{1/.}.bam"
    {% endif %}

{% endif %}


{% if aligner.value == "hisat2" %}

    # Build the Bowtie2 index if needed.
    if [ ! -f "$INDEX.ht2" ]; then
        echo "Building the hisat2 index."
        hisat2-build  ${GENOME} ${INDEX} >> $RUNLOG 2>&1
    fi

    {% if library.value == "SE" %}

        # Run hisat2 in single end mode.
        cat ${FILES} | parallel -j $PROC "hisat2 -x ${INDEX} -U {1} 2>> $RUNLOG | samtools sort > bam/{1/.}.bam"
    {% else %}
        # Run hisat2 in paired end mode.
        cat ${FILES} | parallel -N 2 -j $PROC "hisat2 -x ${INDEX} -1 {1} -2 {2} 2>> $RUNLOG | samtools sort > bam/{1/.}.bam"
    {% endif %}

{% endif %}

# Generate the indices
ls -1 bam/*.bam | parallel samtools index {}

# Generate alignment statistics.
for fname in bam/*.bam; do

    echo "-------- flagstat: $fname -------" >> flagstat.txt
    samtools flagstat ${fname} >> flagstat.txt

    echo "-------- idxstats: $fname -------" >> idxstats.txt
    samtools idxstats ${fname} >> idxstats.txt

done

# Inform the user of the outputs.
echo ""
echo "*** Main Results ***"
echo "Mapping statistics stored in: flagstat.txt"
echo "Alignment counts stored in: idxstats.txt"
echo "********************"
echo ""

# Show a partial output of flagstas.txt
echo ""
echo "**********************"
echo "Partial output for flagstat.txt:"
cat flagstat.txt | head -30
echo "..."
echo "**********************"

# Show a partial output of flagstas.txt
echo ""
echo "**********************"
echo "Partial output for idstats.txt:"
cat idxstats.txt | head -30
echo "..."
echo "**********************"
