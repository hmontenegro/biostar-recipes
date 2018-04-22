"""

"""
import os
import csv

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')




def fetch_unclassified(sample_sheet, source, results, ext='.rep'):

    data = dict()

    stream = csv.DictReader(open(sample_sheet, 'r', newline=''))
    for row in stream:

        result_file = os.path.join(results, row['sample'] + ext)
        # These are the fastq files
        inputs = (os.path.join(source, row['read1']), os.path.join(source, row['read2']))


        data.setdefault(result_file, inputs)

    print(data)

    data = dict()

    1/0

    return data




def main():
    from argparse import ArgumentParser
    import sys

    parser = ArgumentParser()

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

    data = fetch_unclassified(sample_sheet=args.sample_sheet,
                              source=args.input,
                              results=args.results)

    print(data)
    # Check if readID in the sample.rep file is in read1 file or read2 file ( this is where you zip).




if __name__ == '__main__':
    main()

