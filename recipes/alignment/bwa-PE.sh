set -ueo pipefail

# Obtain the parameters.
INPUT={{reads.toc}}
GENOME={{genome.value}}

# Build the BWA index.
INDEX_DIR={{runtime.local_root}}/temp
mkdir -p ${INDEX_DIR}
INDEX=${INDEX_DIR}/{{genome.uid}}

# This directory will store the alignment.
mkdir -p bam

# Build the BWA index it does not already exist.
if [ ! -f "$INDEX.bwt" ]; then
    echo "Building the bwa index."
    bwa index -p ${INDEX} ${GENOME}
else
    echo "Found an existing bwa index."
fi

# Align sequences
echo  "Aligning reads to the genome."
cat ${INPUT} | egrep "fastq|fq" | parallel -N 2 "bwa mem -t 4 ${INDEX} {1} {2} 2>> bwa.log | samtools sort > bam/{1/.}.bam"

# Generate the indices
ls -1 bam/*.bam | parallel samtools index {}

# Generate an alignment report on each.
echo "Compute alignment statistics."

# Reset the file.
echo '' > report.txt

for fname in bam/*.bam; do

    echo "-------- $fname -------" >> report.txt
    samtools flagstat ${fname} >> report.txt

    echo "-------- $fname -------" >> report.txt
    samtools idxstats ${fname} >> report.txt

done

# Show the statistics on the output.
cat report.txt | head -100

