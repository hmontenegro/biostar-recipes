set -ueo pipefail

# Obtain the parameters.
INPUT={{reads.toc}}
GENOME={{genome.value}}
LIBRARY={{library.value}}
ALIGNER={{aligner.value}}
SAMPLING={{sampling.value}}

# Will store the indices here.
# The first time around users must wait until this index is created
# before starting another job that uses this same index.
INDEX_DIR={{runtime.local_root}}/indices
mkdir -p ${INDEX_DIR}
INDEX=${INDEX_DIR}/{{genome.uid}}

# This directory will store the alignment.
mkdir -p bam

# Build the BWA if necessary.
if [ ! -f "$INDEX.bwt" ] && [ "$ALIGNER" == "bwa" ]; then
    echo "Building the bwa index."
    bwa index -p ${INDEX} ${GENOME} >> log.txt 2>&1
fi

# Build the Bowtie2 index if necessary.
if [ ! -f "$INDEX.bt2" ] && [ "$ALIGNER" == "bowtie2" ]; then
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
echo "Compute alignment statistics."

# Reset the file.
echo '' > report.txt

for fname in bam/*.bam; do

    echo "-------- samtools flagstat: $fname -------" >> report.txt
    samtools flagstat ${fname} >> report.txt

    echo "-------- samtools idxstats: $fname -------" >> report.txt
    samtools idxstats ${fname} >> report.txt

done

# Show the statistics on the output.
cat report.txt | head -100

