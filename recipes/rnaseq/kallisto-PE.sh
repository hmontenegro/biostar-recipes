set -ueo pipefail

# Obtain the parameters.
INPUT={{reads.toc}}
TRANSCRIPTS={{transcripts.value}}

# Directory with Kallisto results.
mkdir -p results

# Build the Kallisto index.
INDEX_DIR={{runtime.local_root}}/temp
mkdir -p ${INDEX_DIR}
INDEX=${INDEX_DIR}/{{transcripts.uid}}


# Build the Kallisto index it does not already exist.
if [ ! -f "$INDEX.idx" ]; then
    echo "Building the kallisto index."
    kallisto index -i ${INDEX} ${TRANSCRIPTS}
else
    echo "Found an existing kallisto index."
fi

# Calculate abundances using kallisto.
cat ${INPUT}| egrep "fastq|fq" | sort | parallel -N 2 -j 4 kallisto quant -o results/{1/.}.out  -i ${INDEX} {1} {2}

# Create a combined count table for all samples.

# Get all abundance files.
find results -name "abundance.tsv"  | cut -d "/" -f 2 | parallel cp results/{}/abundance.tsv {}_abundance.txt

# Add sample names to count column.
ls *abundance.txt | parallel sed -i  -e 's/est_counts/{.}/' -e 's/.out_abundance//' {}

# Combine all abundance files.
paste *abundance.txt >all.txt

# Extract count column for all samples.
cat all.txt | awk '{OFS="\t"; for(i=4; i<=NF; i+=5) printf("%s\t",$i); print ""}' >vals.txt

# Extract transcript ids.
cat all.txt | cut -f 1 >ids.txt

# Merge transcript ids with counts to create count table.
paste ids.txt vals.txt > results/counts.txt

# Remove intermediate files.
rm -f *abundance.txt all.txt vals.txt ids.txt *.bak

