# Stop on any error.
set -uxe

# Assembly scaffolds.
SEQUENCE={{sequence.value}}

# Species to be used as training set.
SPECIES={{species.value}}

# No. of processors to be used.
NPROC={{processors.value}}

# Augustus results directory.
AUGUSTUS={{runtime.work_dir}}/augustus

# Create AUGUSTUS directory.
mkdir -p $AUGUSTUS

# Augustus result files.
GENES=${AUGUSTUS}/annotations.gtf
PROTEINS=${AUGUSTUS}/proteins.fa

# Run the Augustus gene predictor in parallel on each contig.
# Use the verbatim  block when a special character is in conflict with
# the Django Templating Language
mkdir -p tmp

{% verbatim %}
cat $SEQUENCE | parallel --j $NPROC --blocksize 5M --recstart '>' --pipe "cat > tmp/{%} && augustus --species=$SPECIES tmp/{%}" > $GENES
{% endverbatim %}
rm -rf tmp

# Make the proteins.fa file from the Augustus predicted GTF file.
getAnnoFasta_mod.pl $GENES

# Change the name of the output.
mv ${AUGUSTUS}/genes.aa ${PROTEINS}



