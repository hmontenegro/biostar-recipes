set -ue

# Obtain the parameters.
INPUT={{reads.toc}}
GENOME={{genome.value}}
LIBRARY={{library.value}}
CUTOFF={{cutoff.value}}

# This is the location of the index directory.
INDEX_DIR=indices
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

# Build the BWA index if needed.
if [ ! -f "$INDEX.bwt" ]; then
    echo "Building the bwa index."
    bwa index -p ${INDEX} ${GENOME} >> $RUNLOG 2>&1
fi

# Run quality trimming.
cat ${FILES} |  parallel -N 2 -j 1 bbduk.sh -Xmx4g maq=20 qtrim=r trimq=15 in1={1} in2={2} out1=${READS}/corrected_{1/} out2=${READS}/corrected_{2/}  overwrite=t

# Merge paired end reads.
echo "Merge paired end reads."
ls -1 ${READS}/corrected_* | parallel -N 2 -j 1 bbmerge.sh ultrastrict=t trimq=20 minoverlap=150 in1={1} in2={2} out=${READS}/merged_{1/} 2>>$RUNLOG

# Run bwa in single end on the merged data.
ls -1 ${READS}/merged_corrected_* | parallel -j $PROC "bwa mem ${INDEX} {1} 2>> $RUNLOG | samtools view -h -q 1 | samtools sort > bam/{1/.}.bam"

# Read stats after merging
echo "Read statistics after merging."
seqkit stat ${READS}/*.fq.gz

# Generate the indices
ls -1 bam/*.bam | parallel samtools index {}
echo "BAM files created in the 'bam' directory."

# Create a results directory
mkdir -p results

# Generate alignment statistics.
ls -1 bam/*bam | parallel "samtools flagstat {} > results/{/}.flagstat.txt"
echo "Flagstats created in the 'results' directory."

ls -1 bam/*bam | parallel "samtools idxstats {} > results/{/}.idxstats.txt"
echo "Index stats created in the 'results' directory."

# Create a combined report of all index stats.
echo "Index statistics in the 'idxstats.txt file"
python -m recipes.code.combine_samtools_idxstats --cutoff $CUTOFF results/*idxstats*t | column -t -s , > idxstats.txt
