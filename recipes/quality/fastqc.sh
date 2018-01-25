set -ue

# The table of contents contains all files in the collection.
INPUTS={{reads.toc}}

# The directory that contains the reports.
mkdir -p results

# Add extra flags if the group parameter was selected.
if [ "{{group.value}}" == "nogroup" ]; then
    FLAGS="--nogroup"
else
    FLAGS=""
fi

# Run fastqc on selected files.
cat ${INPUTS} | egrep 'fastq|fq' | parallel fastqc -q ${FLAGS} {} -o results
