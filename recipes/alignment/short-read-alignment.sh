set -ueo pipefail

# Obtain the parameters.
INPUT={{reads.toc}}
GENOME={{genome.value}}
LIBRARY={{library.value}}
ALIGNER={{aligner.value}}
FRACTION={{fraction.value}}

# Will store the indices here.
# The first time around users must wait until this index is created
# before starting another job that uses this same index.
INDEX_DIR={{runtime.local_root}}/indices
mkdir -p ${INDEX_DIR}
INDEX=${INDEX_DIR}/{{genome.uid}}

# This directory will store the alignment.
mkdir -p bam

# Sub-sample data if necessary.
if [ ${FRACTION} != "1" ]; then
    # Directory with sub-sampled data.
    READS=reads
    mkdir -p $READS

    # A seed of sampling
    SEED=$((1 + RANDOM % 1000))

    # Generate a random sample of each input file.
    echo "Sampling fraction=$FRACTION of data with random seed=$SEED"
    cat ${INPUT} | egrep "fastq|fq" | parallel "seqtk sample -2 -s $SEED {} $FRACTION >$READS/{/.}.fq"

    # The name of the new table of contents.
    INPUT=$READS/toc.txt

    # Create table of contents with sub-sampled data.
    ls -1 $READS/*.fq > ${INPUT}

fi

# Build the BWA index if needed.
if [ "$ALIGNER" == "bwa" ] && [ ! -f "$INDEX.bwt" ]; then
    echo "Building the bwa index."
    bwa index -p ${INDEX} ${GENOME} >> log.txt 2>&1
fi

# Build the Bowtie2 index if needed.
if [ "$ALIGNER" == "bowtie2" ] && [ ! -f "$INDEX.bt2" ]; then
    echo "Building the bowtie2 index."
    bowtie2-build  ${GENOME} ${INDEX} >> log.txt 2>&1
fi

# BWA single end mode.
if [ "$ALIGNER" == "bwa" ] && [ "$LIBRARY" == "SE" ]; then
    echo  "Aligning single end reads to the genome using bwa."
    cat ${INPUT} | sort | egrep "fastq|fq" | parallel "bwa mem -t 4 ${INDEX} {1} 2>> log.txt | samtools sort > bam/{1/.}.bam"
fi

# BWA paired end mode.
if [ "$ALIGNER" == "bwa" ] && [ "$LIBRARY" == "PE" ]; then
    echo  "Aligning paired end reads to the genome using bwa."
    cat ${INPUT} | sort | egrep "fastq|fq" | parallel -N 2 "bwa mem -t 4 ${INDEX} {1} {2} 2>> log.txt | samtools sort > bam/{1/.}.bam"
fi

# Bowtie2 single end mode.
if [ "$ALIGNER" == "bowtie2" ] && [ "$LIBRARY" == "SE" ]; then
    echo  "Aligning single end reads to the genome using bowtie2."
    cat ${INPUT} | sort | egrep "fastq|fq" | parallel "bowtie2 -x ${INDEX} -U {1} 2>> log.txt | samtools sort > bam/{1/.}.bam"
fi

# Bowtie2 paired end mode.
if [ "$ALIGNER" == "bowtie2" ] && [ "$LIBRARY" == "PE" ]; then
    echo  "Aligning paired end reads to the genome using bowtie2."
    cat ${INPUT} | sort | egrep "fastq|fq" | parallel "bowtie2 -x ${INDEX} -1 {1} -2 {2} 2>> log.txt | samtools sort > bam/{1/.}.bam"
fi

# Generate the indices
ls -1 bam/*.bam | parallel samtools index {}

# Generate an alignment report on each.
echo "Computing alignment reports."

# Reset the file.
echo '' > report.txt

for fname in bam/*.bam; do

    echo "-------- flagstat: $fname -------" >> flagstat.txt
    samtools flagstat ${fname} >> flagstat.txt

    echo "-------- idxstats: $fname -------" >> idxstats.txt
    samtools idxstats ${fname} | head -30 >> idxstats.txt

done

# Show the statistics on the output.
cat flagstat.txt | head -100

