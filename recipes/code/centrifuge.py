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
    
    $   python centrifuge.py --files=data/*.txt 
                   
                                                  (fileA, fileB, ...)
        rank-code |  tax-id  |      name        |   percentage-covered  |  reads covered |   reads nassigned
                                                                                                  
        U               0      unclassified        2.00,7.40              1,1              1,7       
        S          314295      Hominoidea          0.01                   5                0
        S          314293      Simiiformes         0.02,.03,.035          1,4,4            1,4,8
    
"""


# Group using the third column ( rank code ).
GROUP_WITH = dict(rank=2)



def generate_header(files):

    header = ["name", "taxID", 'taxRank'] + [os.path.basename(p) for p in files]

    return header

def clean_row(row):
    "Helps clean results of spaces and tabs."

    # Mutates the list
    for i,v in enumerate(row):
        row[i] = v.strip()
    return row


def summarize_group(rank_group):
    "Summarize a rank group by flattening."

    # Take the index of file in header ( exclude the three name, per, rank)

    # dict of index to header name
    summary = []
    for name in rank_group:

        percent = '\t'.join([x[-1] for x in rank_group[name]])
        # Stay constant
        rank = rank_group[name][0][2]
        taxid = rank_group[name][0][1]

        # Match rows to header
        summary.append(f"{name}\t{taxid}\t{rank}\t{percent}")

    return summary


def parse_file(fname, store={}):
    "Parse file and group its contents into a store."

    with open(fname, "r") as infile:
        for row in infile:

            # Inner loop only lasts 6 iterations and only groups the 'rank' column
            for idx, item in enumerate(row.split("\t")):
                if idx == GROUP_WITH['rank']:

                    store.setdefault(item.strip(), []).append(clean_row(row.split("\t")))


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
            if row[0] != 'name':
                name = row[0].strip()
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
        print('\t'.join(generate_header(files=files)))
        for rank, values in summary.items():
            if values:
                print('\n'.join(values))
    else:
        with open(args.outfile, "w") as outfile:
            outfile.write('\t'.join(generate_header(files=files)))
            for rank, values in summary.items():
                if values:
                    outfile.write('\n'.join(values) + "\n")

if __name__ == '__main__':


    try:
        main()
    except (BrokenPipeError, KeyboardInterrupt) as e:
        sys.stderr.close()