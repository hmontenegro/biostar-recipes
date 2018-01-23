set -ue

# The table of content contains all files in the collection.
INPUTS={{reads.toc}}

# Make the directory that contains the reports.
mkdir -p reports

# Add extra flags if a paramter is checked.
if [ "{{group.value}}" == "nogroup" ]; then
    FLAGS="--nogroup"
else
    FLAGS=""
fi

# Run it on fastq files only.
cat ${INPUTS} | egrep 'fastq|fq' | parallel fastqc -q ${FLAGS} {} -o reports
