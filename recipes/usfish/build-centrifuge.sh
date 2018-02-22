# Build fish specific centrifuge database.

ID=7776

ACC=acc.txt
SEQ=sequences.fa
TAXAMAP=taxamap.txt
TAXDIR=/export/refs/taxonomy
INDEX=/export/refs/centrifuge/foo

mkdir -p index

taxonkit list --ids ${ID} --indent "" | sort > targets.txt

#
taxonkit list --ids 40674 --indent "" | sort > remove.txt


# Remove unwanted
comm -2 -3 targets.txt remove.txt  > selected.txt


# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
# centrifuge-download -o $TAXDIR taxonomy

gzcat nucl_gb.accession2taxid.gz | cut -f 2,3 | grep -w -f selected.txt > table.txt

cat table.txt | cut -f 1 > accession.txt

blastdbcmd -db nt -entry_batch accession.txt -out sequence.fa

centrifuge-build -p 4 --conversion-table table.txt --taxonomy-tree $TAXDIR/nodes.dmp --name-table $TAXDIR/names.dmp sequence.fa $INDEX >/dev/null


