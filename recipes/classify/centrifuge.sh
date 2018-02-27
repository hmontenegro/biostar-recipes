
# The table of contents for the data.
INPUT={{reads.toc}}

REFERENCE={{reference.value}}

# The library type.
LIBRARY={{library.value}}

# The sorted file list.
FILES=files.txt

# Create the sorted filelist.
cat ${INPUT}| egrep "fastq|fq" | sort > $FILES

# Select the database to classify against.
if [ ${REFERENCE} == "BAVH" ]; then
    INDEX=/export/refs/centrifuge/p_compressed+h+v
else
    # Custom database for taxid 7776
    INDEX=/export/refs/centrifuge/7776
fi

# Create the reports file.
mkdir -p results

# Choose the right classifier depending on the library setting.
if [ ${LIBRARY} == "PE" ]; then
    # Paired end classification.
    cat ${FILES} | parallel -N 2  -j 1 "centrifuge -x  $INDEX -1 {1} -2 {2} > results/{1/.}.rep"
else
    # Single end classification.
    cat ${FILES} | parallel -j 1 "centrifuge -x  $INDEX -U {} > results/{1/.}.rep"
fi

# Generate an individual report for each sample
for FNAME in results/*.rep; do
    echo "-------- Processing $FNAME -------"
    centrifuge-kreport -x $INDEX $FNAME > $FNAME.txt
done

# Generate a cumulative report as well.
echo "-------- Generating the final report -------"
centrifuge-kreport -x $INDEX results/*.rep > report.txt

echo ""
echo "****************************************"
echo "Individual classification: results"
echo "****************************************"
echo "Cumulative classification in: report.txt"
echo "****************************************"

