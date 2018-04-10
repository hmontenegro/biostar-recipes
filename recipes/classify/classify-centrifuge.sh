set -ue

# The input directory for the data
DDIR=$(dirname {{reads.value}})

# The reference directory to classify against.
REFERENCE={{reference.value}}

# The sum of each row needs to be above this value.
CUTOFF={{cutoff.value}}

# The minimal hit length for classification.
HITLEN={{hitlen.value}}

# The input sample sheet
SHEET={{sheet.value}}

# Output generated while running the tool.
RUNLOG=runlog/runlog.txt

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

# All results go into this folder.
rm -rf results
mkdir -p results

echo "" > $RUNLOG

# How many parallel processes to allow.
PROC=2

# Perform the classification.
cat ${SHEET} | parallel --header : --colsep , -j $PROC  "centrifuge -x  $INDEX -1 $DDIR/{read1} -2 $DDIR/{read2} --min-hitlen $HITLEN -S results/{sample}.rep --report-file  results/{sample}.tsv 2>> $RUNLOG"

set +e
# Generate an individual kraken style reports for each sample
# This will raise an error on no hits, hence we turn of error checking.
ls -1 results/*.rep | parallel -j $PROC "centrifuge-kreport -x $INDEX {} > results/{/.}.txt 2>> $RUNLOG"
set -e

# Generate a combined reformatted.
python -m recipes.code.combine_centrifuge_reports --cutoff $CUTOFF results/*.txt | column -t -s , > classification.txt

