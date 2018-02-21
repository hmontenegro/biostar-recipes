
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
if [ ${REFERENCE} == "HV" ]; then
    # Use the Centrifuge index.
    INDEX=/export/refs/centrifuge/p_compressed+h+v
else
    INDEX=/export/refs/centrifuge/p_compressed+h+v
fi

# Create the reports file.
mkdir -p results

# Choose the right classifier depending on the library setting.
if [ ${LIBRARY} == "PE" ]; then
    # Paired end classification.
    cat ${FILES} | parallel -N 2  -j 4 "centrifuge -x  $INDEX -1 {1} -2 {2} > results/{1/.}.rep"
else
    # Single end classification.
    cat ${FILES} | parallel -j 4 "centrifuge -x  $INDEX -U {} > results/{1/.}.rep"
fi

# Reformat each results as a report.
for FNAME in results/*.rep; do
    echo "-------- Processing $FNAME -------"
    centrifuge-kreport -x $INDEX $FNAME > $FNAME.txt
done

# Rename file to be readable in the browser.
mv centrifuge_report.tsv centrifuge_report.txt
