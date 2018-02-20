
# The table of contents for the data.
INPUT={{reads.toc}}

# The library type.
LIBRARY={{library.value}}

# The sorted file list.
FILES=files.txt

# Create the sorted filelist.
cat ${INPUT}| egrep "fastq|fq" | sort > $FILES

# The Centrifuge index.
INDEX=/export/refs/centrifuge/p_compressed+h+v

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


# Reformat the results as a report.
for FNAME in results/*.rep; do
    echo "-------- Processing $FNAME -------"
    centrifuge-kreport -x $INDEX $FNAME > $FNAME.txt
done


