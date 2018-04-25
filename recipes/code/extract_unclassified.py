import os
import gzip


DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

def decode(stream):
    for elem in stream:
        yield elem.decode("utf8")

def separate(fname, seen=(), outdir=None):
    """
    Takes a fastq stream and a set of already seen/classified reads
    and separate what is not seen/unclassified into a new fastq file.
    """

    if fname.endswith('.gz'):
        instream = decode(gzip.open(fname, 'rb'))
        basename = os.path.splitext(os.path.basename(fname))[0]
    else:
        instream = open(fname, 'rt')
        basename =  os.path.basename(fname)

    outname = "Unclassified_" + basename

    path = os.path.join(outdir, outname)

    outstream = open(path, 'w')

    for rid in instream:
        # Drop leading @ in readid.
        not_classified = rid.split()[0][1:] not in seen

        seq = next(instream)
        tmp = next(instream)
        qual = next(instream)

        if not_classified:
            outstream.write(f"{rid}{seq}{tmp}{qual}")

    return path


def classified_reads(report_files):
    "Get classified reads from report files."

    classified = set()

    for fname in report_files:
        source = open(fname, 'r')
        for line in source:
            classified.add(line.split()[0].strip())

    return classified


def main():
    from argparse import ArgumentParser
    import sys

    parser = ArgumentParser()

    parser.add_argument('files', metavar='FILES', type=str, nargs='+',
                        help='Fastq file(s) used to generate the report file')

    parser.add_argument('--report_files', dest='report', type=str, required=True,
                        nargs='+', help='Generated report file.')

    parser.add_argument('--outdir', dest='outdir', type=str,
                        help='Directory to store the unclassified fastq file.')

    if len(sys.argv) == 1:

        sys.argv.extend([f'{DATA_DIR}/test.fastq.gz', '--report_files', f'{DATA_DIR}/test.rep', f'--outdir={DATA_DIR}'])

    args = parser.parse_args()

    # All classified reads across report files
    seen = classified_reads(report_files=args.report)

    for fname in args.files:

        # Prints the path to unclassified fastq file
        separate(seen=seen, outdir=args.outdir, fname=fname)


if __name__ == '__main__':
    main()

