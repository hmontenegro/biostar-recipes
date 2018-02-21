# Build fish specific centrifuge database.

ID=7776
TAXID=taxid.txt
ACC=acc.txt
SEQ=sequences.fa
TAXAMAP=taxamap.txt
TAXDIR=/export/refs/taxonomy
INDEX=index/centrifuge

mkdir -p index

taxonkit list --ids ${ID} --indent "" > ${TAXID}

# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
# blastdbcmd -entry 'all' -db nt -outfmt '%a,%T' |tr ',' '\t'  >${TAXAMAP}
# centrifuge-download -o $TAXDIR taxonomy


zcat nucl_gb.accession2taxid.gz | \
      csvtk -t grep -f taxid -P ${TAXID} | \
      csvtk -t cut -f gi > ${ACC}

blastdbcmd -db nr -entry all -outfmt "%a\t%T" | \
      csvtk -t grep -f 2 -P ${ACC} | \
      csvtk -t cut -f 1 | \
      blastdbcmd -db nr -entry_batch - -out ${SEQ}

centrifuge-build -p 4 --conversion-table $TAXAMAP --taxonomy-tree $TAXDIR/nodes.dmp --name-table $TAXDIR/names.dmp $SEQ $INDEX >/dev/null
