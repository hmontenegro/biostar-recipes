import os
import sys
import glob


"""

This program is used to process sample report files that are outputted by centrifuge-kreport.
It groups the results according to rank code.

The input files are expected to have six columns:

1. Percentage of reads covered by the clade rooted at this taxon
2. Number of reads covered by the clade rooted at this taxon
3. Number of reads assigned directly to this taxon
4. A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, 
    (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
5. NCBI taxonomy ID
6. indented scientific name


Usage:
   
    $   python centrifuge.py --dir=data/*.txt 
                   
                                                 (fileA, fileB, ...)
        name         |   rank-code   |  tax-id | percentage-covered | reads covered |   reads nassigned
                                                                                                  
        unclassified       U           4            2.00,7.40            1,1            1,7       
        Blochmannia        S           2            0.01                 5              0
        Spirochaeta        S           1            0.02,.03,.035        1,4,4          1,4,8
    
"""


# Group using the third column ( rank code ).
GROUP_WITH = dict(rank=3)

# Header of output
HEADER = "name\trank\ttaxonomy ID\tPercentage of reads\t Number of reads covered\tNumber of reads assigned"


def clean_row(row):
    "Helps clean results of spaces and tabs."

    # Mutates the list
    for i,v in enumerate(row):
        row[i] = v.strip()
    return row


def summarize_group(rank_group):
    "Concatenates information to flatten structure and arranges rows to match header."

    summary = []
    for name in rank_group:

        percent = ','.join([x[0] for x in rank_group[name]])
        ncovered = ','.join([x[1] for x in rank_group[name]])
        nassigned = ','.join([x[2] for x in rank_group[name]])
        # Stay constant
        rank = rank_group[name][0][3]
        taxid = rank_group[name][0][4]

        # Match rows to header
        summary.append(f"{name}\t{rank}\t{taxid}\t{percent}\t{ncovered}\t{nassigned}")

    return summary


def parse_file(fname, store={}):
    "Parse file and group its contents into a store."

    with open(fname, "r") as outfile:
        for row in outfile:

            # Inner loop only lasts 6 iterations and only groups the 'rank' column
            for idx, item in enumerate(row.split("\t")):
                if idx == GROUP_WITH['rank']:
                    val = f"{item.strip()}, rank"
                    store.setdefault(val, []).append(clean_row(row.split("\t")))


def summarize_results(results):
    "Summarize result found in data_dir by grouping them."

    store = dict()
    for item in results:
        if os.path.isfile(item):
            fname = os.path.abspath(item)
            parse_file(fname=fname, store=store)

    for x in store:
        name_store = {}
        for row in store[x]:
            name = row[-1]
            name_store.setdefault(name, []).append(row)
        store[x] = summarize_group(name_store)

    return store


def main():
    from argparse import ArgumentParser

    parse = ArgumentParser()
    parse.add_argument('--files', dest='files', required=True,
                       help="""Directory containing files to be parsed.""",
                       type=str)
    parse.add_argument('--outfile', dest='outfile',
                       help="""Output file to write summary to.""",
                       type=str)
    args = parse.parse_args()

    if len(sys.argv) == 1:
        sys.argv.append('-h')

    files = glob.glob(args.files)
    summary = summarize_results(results=files)

    if not args.outfile:
        print(HEADER)
        for rank in summary:
            print('\n'.join(summary[rank]))
    else:
        with open(args.outfile, "w") as outfile:
            outfile.write(HEADER + "\n")
            for rank in summary:
                outfile.write('\n'.join(summary[rank]) + "\n")

if __name__ == '__main__':


    try:
        main()
    except (BrokenPipeError, KeyboardInterrupt) as e:
        sys.stderr.close()