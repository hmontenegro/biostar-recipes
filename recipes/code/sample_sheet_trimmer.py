"""
Generates the bash commands required to trim data based on a sample sheet.

The output needs to be piped through bash.
"""
import sys, csv, os

FILTER_COMMAND = 'cutadapt --quiet -g ^{fwd} -G ^{rev} --pair-filter both --no-trim --discard-untrimmed {inp_path1} {inp_path2} -o {out_path1} -p {out_path2}'

def generate(args):
    stream = open(args.sheet)
    stream = csv.DictReader(stream, dialect=csv.excel)

    inpdir = args.inpdir
    outdir = args.outdir

    if inpdir:
        inpdir = inpdir.strip("/") + "/"

    if outdir:
        outdir = outdir.strip("/") + "/"

    for row in stream:
        #print (row)
        sample = row['sample']
        fwd = row['barcode'] + row['fwd_primer']
        rev = row['barcode'] + row['rev_primer']
        read1 = row['read1']
        read2 = row['read2']
        min_len = row['min_length']
        max_len = row['max_length']


        name1, ext1 = read1.split(".", 1)
        name2, ext2 = read2.split(".", 1)

        inp_path1 = f"{inpdir}{read1}"
        inp_path2 = f"{inpdir}{read2}"

        if not os.path.isfile(inp_path1) or not os.path.isfile(inp_path2):
            print(f"# Missing files: {inp_path1} or {inp_path2}")
            continue

        out_path1 = f"{outdir}{sample}_1.{ext1}"
        out_path2 = f"{outdir}{sample}_2.{ext2}"

        row.update(dict(inp_path1=inp_path1, inp_path2=inp_path2, out_path1=out_path1, out_path2=out_path2, fwd=fwd, rev=rev))

        filter_command = FILTER_COMMAND.format(**row)

        print(filter_command)


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('sheet', metavar='FILE', type=str,
                        help='The input sample sheet')

    parser.add_argument('--inpdir', type=str, default='',
                        help='The input directory')

    parser.add_argument('--outdir',  type=str, default='',
                        help='The output directory')


    if len(sys.argv) == 1:
        sys.argv.extend([ '--inpdir', 'data', '--outdir', 'out', 'data/sample-sheet.csv'])

    args = parser.parse_args()



    generate(args)


if __name__ == '__main__':
    main()
