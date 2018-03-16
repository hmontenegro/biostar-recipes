import csv

stream = open("names.txt")
stream = csv.reader(stream, delimiter="\t")

out = [ "accession", "taxid", "sciname", "name" ]
print("\t".join(out))
for row in stream:
    row[-1] = row[-1].title()
    print ("\t".join(row))
