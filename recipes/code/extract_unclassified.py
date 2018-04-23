"""

"""
import os
#import csv

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')





def paired_end(stream1, stream2, seen=(), output_dir=None):

    input_stream = zip(stream1, stream2)
    dirname = output_dir or os.path.dirname(stream1.path)

    path1 = os.path.join(dirname,f'Unclassified_{stream1.name}')
    path2 = os.path.join(dirname,f'Unclassified_{stream2.name}')
    output1 = open(path1, 'w')
    output2 = open(path2, 'w')

    for id1, id2, in input_stream:

        seq1, seq2 = next(input_stream)
        tmp1, tmp2 = next(input_stream)
        qual1, qual2 = next(input_stream)

        if id1.strip() not in seen:
            output1.write(f"{id1}{seq1}{tmp1}{qual1}")
        if id2.strip() not in seen:
            output2.write(f"{id2}{seq2}{tmp2}{qual2}")

    return path1, path2


def single_end(stream, seen=(), output_dir=None):

    dirname = output_dir or os.path.dirname(stream.path)
    path = os.path.join(dirname, f'Unclassified_{stream.name}')

    output = open(path, 'w')

    for rid in stream:
        seq = next(stream)
        tmp = next(stream)
        qual = next(stream)

        if id not in seen:
            output.write(f"{rid}{seq}{tmp}{qual}")

    return path


def extract(stream1, report_file, stream2=None):

    seen = set()
    source = open(report_file, 'r')
    for line in source:
        seen.add(line.split()[0])

    if stream2:

        unclassified = [single_end(stream=s, seen=seen) for s in (stream1, stream2)]
    else:
        unclassified = single_end(stream=stream1, seen=seen)


    return unclassified




def main():
    from argparse import ArgumentParser
    import sys

    parser = ArgumentParser()

    parser.add_argument('files', metavar='FILES', type=str, nargs='+',
                        help='Files to be plotted')

    parser.add_argument('--sample_sheet', dest='sample_sheet', type=str, required=True,
                        help='Sample sheet with metabarcoding data.')

    parser.add_argument('--input', dest='input', type=str, required=True,
                        help='Directory with the input reads in Fastq format')

    parser.add_argument('--results', dest='results', type=str, required=True,
                        help='Directory with the classification results')

    if len(sys.argv) == 1:
        sys.argv.extend([f'--sample_sheet={os.path.join(DATA_DIR, "sample-sheet.csv")}',
                         f'--input={DATA_DIR}',
                         f'--results={DATA_DIR}'])

    args = parser.parse_args()


    print(extract())
    # Check if readID in the sample.rep file is in read1 file or read2 file , make new fastq file if not.




if __name__ == '__main__':
    main()

