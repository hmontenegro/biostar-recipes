"""

"""
import os
#import csv

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


def single_end(stream, seen=(), output=None):
    """
    Takes a fastq stream and a set of already seen/classified reads
    and puts what is not seen/unclassified into a fastq file inside of outdir.

    """

    dirname = output or os.path.dirname(stream.path)

    path = os.path.join(dirname, f'Unclassified_{os.path.basename(stream.name)}')

    output = open(path, 'w')

    for rid in stream:
        # Need to further clean line to extract read id.
        check = rid.split()[0].replace('@','')

        seq = next(stream)
        tmp = next(stream)
        qual = next(stream)

        if check not in seen:
            output.write(f"{rid}{seq}{tmp}{qual}")

    return path


def extract(stream1, report_file, stream2=None, output=None):

    seen = set()
    source = open(report_file, 'r')
    for line in source:
        seen.add(line.split()[0].strip())

    if stream2:
        # Handle pair-end reads with list-comp
        unclassified = [single_end(stream=s, seen=seen, output=output) for s in (stream1, stream2)]
    else:
        unclassified = [single_end(stream=stream1, seen=seen, output=output)]

    # Returns list of paths to new fastq files with unclassified reads.
    return unclassified


def main():
    from argparse import ArgumentParser
    import sys

    parser = ArgumentParser()

    parser.add_argument('files', metavar='FILES', type=str, nargs=2,
                        help='Fastq file(s) used to generate the report file')

    parser.add_argument('--report_file', dest='report', type=str, required=True,
                        help='Generated report file.')

    parser.add_argument('--output', dest='output', type=str,
                        help='Directory to store the unclassified fastq file.')

    if len(sys.argv) == 1:
        report_file = os.path.join(DATA_DIR, '1000-MiniXFish.rep')
        files = ['1000-MiniXFish_S2_L001_R1_001.fastq.gz', '1000-MiniXFish_S2_L001_R2_001.fastq.gz']
        files = [os.path.join(DATA_DIR, f) for f in files]

        sys.argv.extend(files + [f'--report_file={report_file}', f'--output={DATA_DIR}'])

    args = parser.parse_args()

    stream1 = open(args.files[0], 'r')
    stream2 = None if len(args.files) <= 1 else open(args.files[1], 'r')

    outputs = extract(stream1=stream1,
                        stream2=stream2,
                        report_file=args.report,
                        output=args.output)
    for o in outputs:
        print(o)


if __name__ == '__main__':
    main()

