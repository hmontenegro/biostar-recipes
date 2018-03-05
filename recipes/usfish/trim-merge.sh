set -ueo pipefail

# Input parameters.

INPUT={{reads.toc}}
SAMPLESHEET={{samplesheet.value}}
QUALITY={{quality.value}}
MIN_LEN={{min_length.value}}
KMER_LEN={{kmer_length.value}}

# Run log to redirect unwanted output.
mkdir -p runlog
RUNLOG=runlog/runlog.txt

# Directory that holds trimmed reads.
mkdir -p trimmed

# Directory that holds merged reads.
mkdir -p merged

# Directory with reads that are filtered for Ns.
mkdir -p filtered

# Get data path.
DATA_DIR=$(dirname $(cat $INPUT| head -1 ))

# Trim forward and reverse primer and then trim by quality.
echo "Removing primers and trimming by quality."

cat $SAMPLESHEET | parallel -j 4  --header : --colsep '\t' bbduk.sh in1=${DATA_DIR}/{read1} in2=${DATA_DIR}/{read2} \
out1=trimmed/{sample}_R1.fq.gz out2=trimmed/{sample}_R2.fq.gz literal={fwd_primer},{rev_primer} ktrim=l k=$KMER_LEN \
hdist=1 qtrim=r trimq=$QUALITY minlength=$MIN_LEN overwrite=true 2>>$RUNLOG

# Merge trimmed reads.
echo "Merging  paired end reads."
ls trimmed/*.fq.gz | sed 's/_R1.fq.gz//g' | parallel -N 2 -j 1 bbmerge.sh in1={1}_R1.fq.gz in2={2} out=merged/{1/}_merged.fq.gz 2>>$RUNLOG

# Remove reads with Ns.
echo " Filtering reads with Ns."
ls merged/*.fq.gz | sed 's/_merged.fq.gz//g' | parallel -j 1 bbduk.sh in={}_merged.fq.gz out=filtered/{/}_filtered.fq.gz maxns=0 2>>$RUNLOG

#
# --------------------------
# Generate a stats table.
# ---------------------------
#
echo "Generating a stats table."

# File with basic stats.
TABLE=stats.txt

echo -e "Sample\tTotal(read1)\tTrimmed\tMerged\tFiltered" >$TABLE

sed 1d $SAMPLESHEET | while IFS='\t' read -r line
do
        sample=$(echo $line| cut -d " " -f 1)
        R1=$(echo $line| cut -d " " -f 7)
        total=$(bioawk -c fastx '{print $name}' $DATA_DIR/$R1 |wc -l)
        trimmed=$(bioawk -c fastx '{print $name}' trimmed/${sample}_R1.fq.gz |wc -l)
        merged=$(bioawk -c fastx '{print $name}' merged/${sample}_merged.fq.gz |wc -l)
        filtered=$(bioawk -c fastx '{print $name}' filtered/${sample}_filtered.fq.gz |wc -l)
        echo -e "$sample\t$total\t$trimmed\t$merged\t$filtered" >>$TABLE
done

