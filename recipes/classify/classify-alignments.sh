set -ue

# Obtain the parameters.
INPUT={{reads.toc}}
GENOME={{genome.value}}
LIBRARY={{library.value}}
CUTOFF={{cutoff.value}}

# This directory holds all intermediary results.
mkdir -p results

# This is the location of the bwa index.
mkdir -p results/index
INDEX=results/index/genome

# Run log to redirect unwanted output.
RUNLOG=runlog/runlog.txt

# Wipe the runlog in case the job is rerun.
echo "" >$RUNLOG

# Build the bwa index."
bwa index -p ${INDEX} ${GENOME} >> $RUNLOG 2>&1

# The name of the input files.
FILES=results/files.txt

# Number of processes.
PROC=4

# This directory should already exist.
mkdir -p runlog

# Generate the input file names.
cat ${INPUT} | sort | egrep "fastq|fq" > ${FILES}


# Error correct the sequences.
mkdir -p results/corrected
cat ${FILES} | parallel -N 2 -j 1 tadpole.sh in1={1} in2={2} out1=results/corrected/{1/} out2=results/corrected/{2/}  mode=correct k=50 overwrite=t 2>>$RUNLOG

# Trim sequences by quality.
mkdir -p results/trimmed
ls -1 results/corrected/* |  parallel -N 2 -j 1 bbduk.sh -Xmx4g maq=20 qtrim=r trimq=15 in1={1} in2={2} out1=results/trimmed/{1/} out2=results/trimmed/{2/}  overwrite=t 2>>$RUNLOG

# Merge trimmed reads.
mkdir -p results/merged
ls -1 results/trimmed/* | parallel -N 2 -j 1 bbmerge.sh ultrastrict=t trimq=20 minoverlap=150 in1={1} in2={2} out=results/merged/{1/} 2>>$RUNLOG

# Read stats after merging
echo "--- Corrected --- "
seqkit stat results/corrected/*
echo "--- Trimmed --- "
seqkit stat results/trimmed/*
echo "--- Merged --"
seqkit stat results/merged/*

# Remove the unused results
# rm -rf results/corrected results/trimmed

# Run bwa on the merged results.
mkdir -p results/bam
ls -1 results/merged/* | parallel -j $PROC "bwa mem ${INDEX} {1} 2>> $RUNLOG | samtools view -h -q 1 | samtools sort > results/bam/{1/.}.bam"
ls -1 results/bam/*.bam | parallel samtools index {}

# Create a results directory
mkdir -p results/counts

# Generate alignment statistics.
ls -1 results/bam/*bam | parallel "samtools flagstat {} > results/counts/{/}.flagstat.txt"
ls -1 results/bam/*bam | parallel "samtools idxstats {} > results/counts/{/}.idxstats.txt"

# Create a combined report of all index stats.
echo "Final results in the 'idxstats.txt file"
python -m recipes.code.combine_samtools_idxstats --cutoff $CUTOFF results/counts/*idxstats*t | column -t -s , > idxstats.txt
