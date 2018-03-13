import csv

def flatten():

    fname = 'species_list.txt'

    stream = csv.reader(open(fname, 'rt'), delimiter='\t')
    next(stream)

    for row in stream:
        common, sci = row[:2]
        rest = list(filter(None, row[2:]))
        for elem in rest:
            out = [common, sci, elem]
            print("\t".join(out))


def extract():
    fname = 'nucl_gb.accession2taxid'

    fname = 'species.txt'
    stream = csv.reader(open(fname, 'rt'), delimiter='\t')
    next(stream)
    keep = set( row[2] for row in stream)

    #print (keep)

    fname = 'nucl_gb.accession2taxid'
    stream = csv.reader(open(fname, 'rt'), delimiter='\t')
    row = next(stream)
    for row in stream:
        acc = row[0]
        if acc in keep:
            out = [row[0], row[2] ]
            print ("\t".join(out))
            #break

if __name__ == '__main__':
    #flatten()

    extract()
