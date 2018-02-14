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

# This directory should already exist.
mkdir -p runlog

# Generate the input file names.
cat ${INPUT} | sort | egrep "fastq|fq" > ${FILES}

# Sub-sample data if necessary.
if [ ${FRACTION} != "1" ]; then

    # The random seed for sampling.
    SEED=$((1 + RANDOM % 1000))

    # Generate a random sample of each input file.
    echo "Sampling fraction=$FRACTION of data with random seed=$SEED"
    cat ${FILES} | egrep "fastq|fq" | parallel "seqtk sample -2 -s $SEED {} $FRACTION >$READS/{/.}.fq"

    # Create a new table of contents with sub-sampled data.
    ls -1 $READS/*.fq > ${FILES}

fi

# Build the BWA index if needed.
if [ "$ALIGNER" == "bwa" ] && [ ! -f "$INDEX.bwt" ]; then
    echo "Building the bwa index."
    bwa index -p ${INDEX} ${GENOME} >> $RUNLOG 2>&1
fi

# Build the Bowtie2 index if needed.
if [ "$ALIGNER" == "bowtie2" ] && [ ! -f "$INDEX.bt2" ]; then
    echo "Building the bowtie2 index."
    bowtie2-build  ${GENOME} ${INDEX} >> $RUNLOG 2>&1
fi

# Build the Hisat2 index if needed.
if [ "$ALIGNER" == "hisat2" ] && [ ! -f "$INDEX.ht2" ]; then
    echo "Building the hisat2 index."
    hisat2-build  ${GENOME} ${INDEX} >> $RUNLOG 2>&1
fi

# BWA single end mode.
if [ "$ALIGNER" == "bwa" ] && [ "$LIBRARY" == "SE" ]; then
    echo  "Aligning single end reads to the genome using $ALIGNER."
    cat ${FILES} | parallel "bwa mem -t 4 ${INDEX} {1} 2>> $RUNLOG | samtools sort > bam/{1/.}.bam"
fi

# BWA paired end mode.
if [ "$ALIGNER" == "bwa" ] && [ "$LIBRARY" == "PE" ]; then
    echo  "Aligning paired end reads to the genome using $ALIGNER."
    cat ${FILES} | parallel -N 2 "bwa mem -t 4 ${INDEX} {1} {2} 2>> $RUNLOG | samtools sort > bam/{1/.}.bam"
fi

# Bowtie2 single end mode.
if [ "$ALIGNER" == "bowtie2" ] && [ "$LIBRARY" == "SE" ]; then
    echo  "Aligning single end reads to the genome using $ALIGNER."
    cat ${FILES} | parallel "bowtie2 -x ${INDEX} --sensitive-local -U {1} 2>> $RUNLOG | samtools sort > bam/{1/.}.bam"
fi

# Bowtie2 paired end mode.
if [ "$ALIGNER" == "bowtie2" ] && [ "$LIBRARY" == "PE" ]; then

    cat ${FILES} | parallel "bowtie2 -x ${INDEX} --sensitive-local -1 {1} -2 {2} 2>> $RUNLOG | samtools sort > bam/{1/.}.bam"
fi

# Hisat2 single end mode.
if [ "$ALIGNER" == "hisat2" ] && [ "$LIBRARY" == "SE" ]; then
    echo  "Aligning single end reads to the genome using $ALIGNER."
    cat ${FILES} | parallel "hisat2 -x ${INDEX} -U {1} 2>> $RUNLOG | samtools sort > bam/{1/.}.bam"
fi

# Hisat2 paired end mode.
if [ "$ALIGNER" == "hisat2" ] && [ "$LIBRARY" == "PE" ]; then
    echo  "Aligning paired end reads to the genome using $ALIGNER."
    cat ${FILES} | parallel "hisat2 -x ${INDEX} -1 {1} -2 {2} 2>> $RUNLOG | samtools sort > bam/{1/.}.bam"
fi

# Generate the indices
ls -1 bam/*.bam | parallel samtools index {}

# Generate alignment statistics.
for fname in bam/*.bam; do

    echo "-------- flagstat: $fname -------" >> flagstat.txt
    samtools flagstat ${fname} >> flagstat.txt

    echo "-------- idxstats: $fname -------" >> idxstats.txt
    samtools idxstats ${fname} | head -30 >> idxstats.txt

done

# Show the statistics on the output.
cat flagstat.txt | head -100

# Generate IGV plots for the alignments
# This step is optional and is used to demonstrate
# a utility script that come with the biostar-recipes

# Make a directory for the gnome
mkdir -p refs

# This will be a link to the genome in job directory.
LOCAL=refs/genome.fa

# Link the genome to the current directory
ln -sf $GENOME $LOCAL

# Create an index for the genome
samtools faidx $LOCAL

# This module crawls the current job directory
# and generates an IGV entry for eligible files in it.
python -m recipes.code.igv --genome  $LOCAL --baseurl {{runtime.job_url}} > igv.xml
