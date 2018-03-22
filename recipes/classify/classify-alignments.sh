set -uex

# Obtain the parameters.
INPUT={{reads.toc}}
GENOME={{genome.value}}
LIBRARY={{library.value}}
CUTOFF={{cutoff.value}}
SHEET={{sheet.value}}

MINLEN={{minlen.value}}

rm -rf results

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

# Filter by adapter sequences. This must be performed by the information in a sample sheet.
mkdir -p results/filtered

# Save the trimming commands to a script for reference.
python -m recipes.code.sample_sheet_trimmer --inpdir results/corrected --outdir results/filtered ${SHEET} > results/trimming.sh

# Run the trimming commands
cat results/trimming.sh | bash

# Merge corrected, filtered reads.
mkdir -p results/merged
ls -1 results/filtered/* | parallel -N 2 -j 1 bbmerge.sh ultrastrict=t minoverlap=$MINLEN in1={1} in2={2} out=results/merged/{1/} 2>>$RUNLOG

# Read stats after merging
echo "--- Corrected --- "
seqkit stat results/corrected/*
echo "--- Filtered --- "
seqkit stat results/filtered/*
echo "--- Merged --"
seqkit stat results/merged/*

# Generate the alignments for each data
# rm -rf results/corrected results/trimmed
mkdir -p results/bam

# These steps are not strictly necessary but are here for troubleshootig for now.
cat ${FILES} | parallel -N 2 -j $PROC "bwa mem ${INDEX} {1} {2} 2>> $RUNLOG | samtools sort > results/bam/original-{1/.}.bam"
ls -1 results/corrected/* | parallel -N 2 -j $PROC "bwa mem ${INDEX} {1} {2} 2>> $RUNLOG | samtools sort > results/bam/corrected-{1/.}.bam"
ls -1 results/filtered/* | parallel -N 2 -j $PROC "bwa mem ${INDEX} {1} {2} 2>> $RUNLOG | samtools sort > results/bam/filtered-{1/.}.bam"

# Generate the merged alignment.
ls -1 results/merged/* | parallel -j $PROC "bwa mem ${INDEX} {} 2>> $RUNLOG | samtools view -h -q 1 -F 2304 | python -m recipes.code.bamfilter --minlen $MINLEN | samtools sort > results/bam/merged-{/.}.bam"

# Generate the indices
ls -1 results/bam/*.bam | parallel samtools index {}

# Create a results directory
mkdir -p results/counts

# Generate alignment statistics.
ls -1 results/bam/merged-*.bam | parallel "samtools flagstat {} > results/counts/{/}.flagstat.txt"
ls -1 results/bam/merged-*.bam | parallel "samtools idxstats {} > results/counts/{/}.idxstats.txt"

# Create a combined report of all index stats.
echo "Final results in the 'idxstats.txt file"
python -m recipes.code.combine_samtools_idxstats --cutoff $CUTOFF results/counts/*idxstats*t | column -t -s , > idxstats.txt
