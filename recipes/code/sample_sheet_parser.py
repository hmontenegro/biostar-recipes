"""
Reads a sample sheet and generates file names and primers as a tabular text file.

Output format has 6 tab delimited columns

read1 read2 primer1 primer2 sample

"""
import sys, csv, os

def parse(args):
    stream = open(args.sheet)
    stream = csv.DictReader(stream, dialect=csv.excel)

    for row in stream:

        sample = row['sample']

        fwd = row['barcode'] + row['fwd_primer']
        rev = row['barcode'] + row['rev_primer']

        read1 = row['read1']
        read2 = row['read2']

        line = [ read1, read2, fwd, rev, sample]

        print ("\t".join(line))


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('sheet', metavar='FILE', type=str,
                        help='The input sample sheet')

    if len(sys.argv) == 1:
        sys.argv.extend(['data/sample-sheet.csv'])

    args = parser.parse_args()

    parse(args)


if __name__ == '__main__':
    main()
