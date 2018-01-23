set -ueo pipefail

# Obtain the parameters.
INPUT={{reads.value}}
GENOME={{genome.value}}

# Build the BWA index
mkdir -p index
INDEX=index/{{genome.uid}}

# This directory will store the alignment.
mkdir -p bam
BAM=bam/alignment.bam

# Build the BWA index.
echo "Build the bwa index."
bwa index -p ${INDEX} ${GENOME}

# Align sequences
echo  "Aligning reads to the genome."
bwa mem -t 4 ${INDEX} ${INPUT}  | samtools view -b |samtools sort >${BAM}
samtools index ${BAM}

echo "Compute alignment statistics."
samtools flagstat ${BAM} > flagstat.txt
samtools idxstats ${BAM} > idxstats.txt

# Show the statistics on the output.
echo "-------- Flag Stats -------------"
cat flagstat.txt

echo "-------- Index Stats ------------"
cat idxstats.txt | head -100
