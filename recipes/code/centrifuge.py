import os
import sys


"""

This program is used to process files that are outputted by centrifuge-kreport.

The files are expected to have six columns:

1. Percentage of reads covered by the clade rooted at this taxon
2. Number of reads covered by the clade rooted at this taxon
3. Number of reads assigned directly to this taxon
4. A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, 
    (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
5. NCBI taxonomy ID
6. indented scientific name


Usage:

    python centrifuge.py --dir data --group_by=rank --output report.txt    # Writes results to a file
    python centrifuge.py --dir data --group_by=rank                        # Outputs to screen

Sample Output -- grouped by the rank code ( column 4 ):

    2.00    1   1   U   4   unclassified
    7.40    1   7   U   0   unclassified
    0.70    4   3   S   14 Spirochaeta thermophila
    0.05    3   7   S   154 Spirochaeta 
    1.01    1   4   G   21   Candidatus Blochmannia
    0.01    5   0   G   2   Blochmannia
    

"""


# Options provided when it comes to grouping results
GROUPING_CHOICES = (

        ("percent", 0) ,
        ("ncovered", 1),
        ("nassigned", 2),
        ("rank", 3),
        ("taxid", 4),
        ("name", 5)
)

def clean_row(row):
    "Helps clean results of spaces and tabs."

    # Mutates the list
    for i,v in enumerate(row):
        row[i] = v.strip()
    return row


def parse_file(rep_file, store={}):
    "Parse a centrifuge-kreport output file and populate dict with contents."

    idx_to_name = {y:x for x,y in GROUPING_CHOICES}
    with open(rep_file, "r") as outfile:
        for row in outfile:
            # Inner loop only lasts 6 iterations
            for idx, item in enumerate(row.split("\t")):

                # Store a string that contains an items group.
                val = f"{item.strip()}, {idx_to_name[idx]}"
                store.setdefault(val, []).append(clean_row(row.split("\t")))


def summarize_results(data_dir, group_by="rank"):
    "Summarize result found in data_dir by grouping them."

    store = dict()
    for item in os.scandir(data_dir):
        if item.is_file():
            fname = os.path.abspath(item.path)
            parse_file(rep_file=fname, store=store)

    summary = []
    for x in store:
        if group_by in x:
            summary.extend(store[x])
    return summary


def main():
    from argparse import ArgumentParser

    parse = ArgumentParser()
    parse.add_argument('--dir', dest='indir', required=True,
                       help="""Directory containing .rep/.rep.txt files to be parsed.""",
                       type=str)
    parse.add_argument('--group_by', dest='group_by', default="rank",
                       help="""Group resulting report in specific manner.""",
                       type=str)
    parse.add_argument('--outfile', dest='outfile',
                       help="""Output file to write summary to.""",
                       type=str)
    args = parse.parse_args()
    summary = []

    if len(sys.argv) == 1:
        sys.argv.append('-h')

    if os.path.isdir(args.indir):
        summary = summarize_results(data_dir=args.indir, group_by=args.group_by)
    else:
        parse.print_help()
        print('--dir needs to be a valid directory.')
        exit()

    if not args.outfile:
        for row in summary:
            print('\t'.join(row))
    else:
        with open(args.outfile, "w") as outfile:
            for row in summary:
                outfile.write('\t'.join(row) + "\n")

if __name__ == '__main__':


    try:
        main()
    except (BrokenPipeError, KeyboardInterrupt) as e:
        sys.stderr.close()