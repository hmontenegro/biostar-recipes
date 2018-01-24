set -ueo pipefail

# Obtain the parameters.
INPUT={{reads.value}}
GENOME={{genome.value}}

# Build the BWA index.
INDEX_DIR={{runtime.local_root}}/temp
mkdir -p ${INDEX_DIR}
INDEX=${INDEX_DIR}/{{genome.uid}}

# This directory will store the alignment.
mkdir -p bam
BAM=bam/alignment.bam

# Build the BWA index it does not already exist.
if [ ! -f "$INDEX.bwt" ]; then
    echo "Building the bwa index."
    bwa index -p ${INDEX} ${GENOME}
else
    echo "Found an existing bwa index."
fi

# Align sequences
echo  "Aligning reads to the genome."
bwa mem -t 4 ${INDEX} ${INPUT}  | samtools sort > ${BAM}
samtools index ${BAM}

# Produce an alignment report.
echo "Compute alignment statistics."
echo "-------- $fname -------" > report.txt
samtools flagstat ${BAM} >> report.txt
echo "-------- $fname -------" >> report.txt
samtools idxstats ${BAM} >> report.txt

# Show the statistics on the output.
cat report.txt | head -100
