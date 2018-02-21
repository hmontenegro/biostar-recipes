# Input parameters.

INPUT={{reads.value}}
SAMPLESHEET={{samplesheet.value}}
QUALTIY={{quality.value}}
MIN_LEN={{min_length.value}}
KMER_LEN={{kmer_length.value}}

# Run log to redirect unwanted output.
mkdir -p runlog
RUNLOG=runlog/runlog.txt

# Directory that holds trimmed reads.
mkdir -p trimmed

# Directory that holds merged reads.
mkdir -p merged

# Directory with reads that are filtered for Ns
mkdir -p filtered

# trim forward and reverse primer and then trim by quality.

cat $SAMPLESHEET | parallel -j 4  --header : --colsep '\t' bbduk.sh in1={read1} in2={read2}  out1=trimmed/{sample}_R1.fq.gz out2=trimmed/{sample}_R2.fq.gz \
literal={fwd_primer},{rev_primer} ktrim=l k=$KMER_LEN hdist=1 qtrim=r trimq=$QUALITY minlength=$MIN_LEN overwrite=true 2>>$RUNLOG


# Merge trimmed reads.
echo "Merging  paired end reads."
ls trimmed/*.fq.gz | sed 's/_R1.fq.gz//g' | parallel -N 2 -j 4 bbmerge.sh in1={1}_R1.fq.gz in2={2} out=merged/{1/}_merged.fq.gz 2>>$RUNLOG

# Remove reads with Ns.
echo "Filtering reads with Ns."
ls merged/*.fq.gz | sed 's/_merged.fq.gz//g' | parallel -j 4 bbduk.sh in={}_merged.fq.gz out=filtered/{}_filtered.fq.gz maxns=0 2>>$RUNLOG


