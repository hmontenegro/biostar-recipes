
# The table of contents for the data.
INPUT={{reads.toc}}

# The Centrifuge index.
INDEX=/export/refs/centrifuge/p_compressed+h+v

# Create the reports file.
mkdir -p results

# Perform a classification on each sequencing file.
cat ${INPUT} | egrep "fastq|fq" | parallel -j 1 "centrifuge -x  $INDEX {} > results/{1/.}.rep"

# Reformat the results as a report.
for FNAME in results/*.rep; do
    echo "-------- Processing $FNAME -------"
    centrifuge-kreport -x $INDEX $FNAME > $FNAME.txt
done


