import os
import sys
import gzip


DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


def join(*paths):
    return os.path.abspath(os.path.join(*paths))


def separate(fname, seen=(), outdir=None):
    """
    Takes a fastq stream and a set of already seen/classified reads
    and separate what is not seen/unclassified into a new fastq file.
    """

    is_gz = fname.endswith('.gz')
    create_stream = lambda p,m: gzip.open(p,m) if is_gz else open(p,m)

    instream = create_stream(fname, 'rt')

    basename = os.path.splitext(os.path.basename(fname))[0]
    basename = basename if is_gz else os.path.basename(fname)
    outname = "Unclassified_" + basename

    path = join(str(outdir), outname)
    outstream = sys.stdout if not outdir else create_stream(path, 'w')

    for rid in instream:
        # Drop leading @ in readid.
        not_classified = rid.split()[0][1:] not in seen

        seq = next(instream)
        tmp = next(instream)
        qual = next(instream)

        if not_classified:
            outstream.write(f"{rid}{seq}{tmp}{qual}")

    return path


def get_classified_reads(report_files):
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

        sys.argv.extend([join(DATA_DIR,"test.fastq.gz"), '--report_files', join(DATA_DIR,'test.rep')])

    args = parser.parse_args()

    # Classified reads across report files
    seen = get_classified_reads(report_files=args.report)

    for fname in args.files:
        separate(seen=seen, outdir=args.outdir, fname=fname)


if __name__ == '__main__':
    main()

