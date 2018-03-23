set -uex

# The data names get sorted lexicographically.
# Hence an explicit listing of data is usually required.
# This will be achieved via a sample file.

# We only need the directoy where the data resides.

DDIR=$(dirname {{reads.value}})

GENOME={{genome.value}}

CUTOFF={{cutoff.value}}

SHEET={{sheet.value}}

MINLEN={{minlen.value}}

# This directory will hold all intermediary results.
rm -rf results

# This directory holds all intermediary results.
mkdir -p results

# This is the location of the bwa index.
mkdir -p results/index
INDEX=results/index/genome

# Number of processes.
PROC=4

# Run log to redirect unwanted output.
RUNLOG=runlog/runlog.txt

# Wipe the runlog in case the job is rerun.
echo "" >$RUNLOG

# Build the bwa index."
bwa index -p ${INDEX} ${GENOME} >> $RUNLOG 2>&1

# The name of the input files.
FILES=results/files.txt

# Generates the six column tab delimited file that contains
# file1 file2 primer1 primer2 output1 output2 output3
python -m recipes.code.sample_sheet_parser $SHEET > $FILES

# This directory should already exist.
mkdir -p runlog

# Error correct the sequences.
CDIR=results/corrected
mkdir -p $CDIR
cat ${FILES} | parallel --colsep '\t' -j 2 tadpole.sh in1=$DDIR/{1} in2=$DDIR/{2} out1=$CDIR/{1} out2=$CDIR/{2}  mode=correct k=50 overwrite=t 2>>$RUNLOG

# Use cutadapt to filter out those read that DO NOT contain the primers
FDIR=results/filtered
mkdir -p $FDIR
cat ${FILES} | parallel --colsep '\t' -j 2 cutadapt --quiet -g ^{3} -G ^{4} --pair-filter both --no-trim --discard-untrimmed $CDIR/{1} $CDIR/{2} -o $FDIR/{1} -p $FDIR/{2}

# Merge corrected, filtered reads. The 7th column is the sample name plus extension.
MDIR=results/merged
mkdir -p $MDIR
cat ${FILES} | parallel --colsep '\t' -j 2 bbmerge.sh ultrastrict=t minoverlap=$MINLEN in1=$FDIR/{1} in2=$FDIR/{2} out=$MDIR/{5}.fq 2>>$RUNLOG

# Read stats after merging
echo "--- Corrected --- "
seqkit stat $CDIR/*
echo "--- Filtered --- "
seqkit stat $FDIR/*
echo "--- Merged --"
seqkit stat $MDIR/*

# Generate the alignments for each data
# rm -rf $CDIR $FDIR

BDIR=results/bam
mkdir -p $BDIR

# These steps are optional and needed only when investigating/debugging.
cat ${FILES} | parallel --colsep '\t' -j 2 "bwa mem ${INDEX} $DDIR/{1} $DDIR/{2} 2>> $RUNLOG | samtools sort > $BDIR/original-{5}.bam"
cat ${FILES} | parallel --colsep '\t' -j 2 "bwa mem ${INDEX} $CDIR/{1} $CDIR/{2} 2>> $RUNLOG | samtools sort > $BDIR/corrected-{5}.bam"
cat ${FILES} | parallel --colsep '\t' -j 2 "bwa mem ${INDEX} $FDIR/{1} $FDIR/{2} 2>> $RUNLOG | samtools sort > $BDIR/filtered-{5}.bam"

# Generate the merged alignment that will be used.
cat ${FILES} | parallel --colsep '\t' -j 2 "bwa mem ${INDEX} $MDIR/{5}.fq 2>> $RUNLOG | samtools view -h -q 1 -F 2304 | python -m recipes.code.bamfilter --minlen $MINLEN | samtools sort > $BDIR/{5}.bam"

# Generate the indices for all BAM files.
ls -1 $BDIR/*.bam | parallel samtools index {}

# Create a results directory
mkdir -p results/counts

# Generate alignment statistics.
cat ${FILES}  | parallel --colsep '\t' "samtools flagstat $BDIR/{5}.bam > results/counts/{5}.flagstat.txt"
cat ${FILES}  | parallel --colsep '\t' "samtools idxstats $BDIR/{5}.bam > results/counts/{5}.idxstats.txt"

# Create a combined report of all index stats.
echo "Final results in the 'idxstats.txt file"
python -m recipes.code.combine_samtools_idxstats --cutoff $CUTOFF results/counts/*idxstats*t | column -t -s , > idxstats.txt
