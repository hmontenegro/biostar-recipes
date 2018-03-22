set -u

# The table of contents for the data.
INPUT={{reads.toc}}

REFERENCE={{reference.value}}

# The library type.
LIBRARY={{library.value}}

# The sum of each row needs to be above this value.
CUTOFF={{cutoff.value}}


# The minimal hit length for classification.
HITLEN={{hitlen.value}}

# The sorted file list.
FILES=runlog/files.txt

# Create the sorted filelist.
cat ${INPUT}| egrep "fastq|fq" | sort > $FILES

# Select the database to classify against.
if [ ${REFERENCE} == "BAVH" ]; then
    INDEX=/export/refs/centrifuge/p_compressed+h+v
fi

if [ ${REFERENCE} == "JAW" ]; then
    # Custom database for taxid 7776
    INDEX=/export/refs/centrifuge/7776
fi

if [ ${REFERENCE} == "FISH" ]; then
    # Custom database for taxid 7776
    INDEX=/export/refs/centrifuge/fishdb
fi

if [ ${REFERENCE} == "NT" ]; then
    # Refseq NT
    INDEX=/export/refs/centrifuge/nt
fi

# Create the reports file.
mkdir -p results

echo "Individual reports saved into the results folder"
# Choose the right classifier depending on the library setting.
if [ ${LIBRARY} == "PE" ]; then
    # Paired end classification.
    cat ${FILES} | parallel -N 2  -j 1 "centrifuge -x  $INDEX -1 {1} -2 {2} --min-hitlen $HITLEN -S results/{1/.}.rep --report-file  results/{1/.}.tsv"
else
    # Single end classification.
    cat ${FILES} | parallel -j 1 "centrifuge -x  $INDEX -U {} --min-hitlen $HITLEN -S results/{1/.}.rep --report-file results/{1/.}.tsv"
fi


# Generate an individual kraken report for each sample
ls -1 results/*.rep | parallel -j 1 "centrifuge-kreport -x $INDEX {} > results/{1/.}.txt"

# Generate a combined reformatted.
echo "Combined report: classification.txt"
echo "Combined plot: classification.png"
python -m recipes.code.combine_centrifuge_reports --cutoff $CUTOFF --plot classification.png results/*.txt | column -t -s , > classification.txt




