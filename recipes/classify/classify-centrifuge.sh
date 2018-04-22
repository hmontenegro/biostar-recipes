set -ue

# The input directory for the data
DDIR=$(dirname {{reads.value}})

# The reference directory to classify against.
REFERENCE={{reference.value}}

# The sum of each row needs to be above this value.
CUTOFF={{cutoff.value}}

# The minimal hit length for classification.
HITLEN={{hitlen.value}}

# The list of files in the directory.
SHEET={{sheet.value}}

# Output generated while running the tool.
RUNLOG=runlog/runlog.txt

# The name of the centrifuge index.
INDEX=index/database

# Create the custom database.
mkdir -p index

# All results from centrifuge go into this folder.
rm -rf results
mkdir -p results

# Wipe clean the runlog.
echo "" > $RUNLOG

# How many parallel processes to allow.
N=2

# Download taxonomy specific files. Run these in  $TAXDIR.
# This operation only needs to be done once for the entire website.
#
# (cd $TAXDIR && wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)
# (cd $TAXDIR && wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz)
# (cd $TAXDIR && gunzip taxdump.tar.gz)
# (cd $TAXDIR && gunzip nucl_gb.accession2taxid.gz)
# Create the conversion table (accession to taxid mapping).
# cat $TAXDIR/nucl_gb.accession2taxid | cut -f 1,3 > $TAXDIR/table.txt
# Untar file
# tar -xvf $TAXDIR/taxdump.tar

# Use the taxonomy specific files to build the custom database.
TAXDIR=/export/refs/taxonomy
TABLE=$TAXDIR/table.txt
NODES=$TAXDIR/nodes.dmp
NAMES=$TAXDIR/names.dmp

# Directory to store classification results
CLASSDIR=classification
rm -rf $CLASSDIR
mkdir -p $CLASSDIR

# Build the index.
centrifuge-build -p $N --conversion-table $TABLE --taxonomy-tree $NODES  --name-table $NAMES  $REFERENCE $INDEX >> $RUNLOG

# Perform the classification.
cat ${SHEET} | parallel --header : --colsep , -j $N  "centrifuge -x  $INDEX -1 $DDIR/{read1} -2 $DDIR/{read2} --min-hitlen $HITLEN -S results/{sample}.rep --report-file  results/{sample}.tsv 2>> $RUNLOG"

set +e
# Generate an individual kraken style reports for each sample
# This will raise an error on no hits, hence we turn of error checking.
ls -1 results/*.rep | parallel -j $N "centrifuge-kreport -x $INDEX {} > results/{/.}.txt 2>> $RUNLOG"
set -e

# Generate a combined reformatted report inside of CLASSDIR.
python -m recipes.code.combine_centrifuge_reports --cutoff $CUTOFF results/*.txt --outdir $CLASSDIR

# Draw the heat maps for each csv report
python -m recipes.code.plotter $CLASSDIR/*.csv --type heat_map

# Generate report for any unclassified reads
python -m recipes.code.unclassified --sample_sheet $SHEET --input_dir $DDIR/.. --results_dir $DDIR/../results

# Draw the rarefaction curves.
python -m recipes.code.rarefaction results/*.rep

