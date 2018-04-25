import os

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')



def single_end(stream, seen=(), output=None):
    """
    Takes a fastq stream and a set of already seen/classified reads
    and puts what is not seen/unclassified into a new fastq file.

    """

    dirname = output or os.path.dirname(stream.name)

    path = os.path.join(dirname, f'Unclassified_{os.path.basename(stream.name)}')

    output = open(path, 'w')

    for rid in stream:
        # Need to further clean line to extract read id.
        not_classified = rid.split()[0].replace('@','') not in seen

        seq = next(stream)
        tmp = next(stream)
        qual = next(stream)

        if not_classified:
            output.write(f"{rid}{seq}{tmp}{qual}")

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

    parser.add_argument('--output', dest='output', type=str,
                        help='Directory to store the unclassified fastq file.')

    if len(sys.argv) == 1:

        sys.argv.extend([f'{DATA_DIR}/*.fastq.gz' ,'--report_files', f'{DATA_DIR}/*.rep', f'--output={DATA_DIR}'])

    args = parser.parse_args()

    # All classified reads across report files
    seen = classified_reads(report_files=args.report)

    for fname in args.files:
        stream = open(fname, 'r')

        # Prints the path to unclassified fastq file
        print(single_end(seen=seen, output=args.output, stream=stream))


if __name__ == '__main__':
    main()

