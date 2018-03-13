"""

This program is used to process .tsv report files that are outputted by centrifuge-kreport.
It groups the results according to rank code and shows the abundance for each file.

The input files are expected to have this header:

    name	taxID	taxRank	genomeSize	numReads	numUniqueReads	abundance

Usage:
    
    $   python centrifuge.py --files=data/*.tsv

    name	                    taxID	taxRank	    SRR1972972_1.tsv	SRR1972972_2.tsv    ... 
    Azorhizobium caulinodans	7	    species	    0.0	                0.0                 
    Cellulomonas gilvus	        11	    species	    -                   0.0
    Pelobacter carbinolicus	    19	    species	    0.0                 -
    ...
                   
"""

import os
import sys
import glob



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


def summarize_group(rank_group, file_columns):
    "Summarize a rank group by flattening."
    
    summary = []
    col_to_idx = {x:y for y,x in enumerate(file_columns)}

    for name in rank_group:

        # Start off with empty values and substitute
        # every column that has a value
        percents = ["-" for _ in range(len(file_columns))]
        current = [x[-1] for x in rank_group[name]]
        for p in current:
            percent, col = p.split('\t')[0], p.split('\t')[1]
            percents[col_to_idx[col]] = percent

        percents = '\t'.join(percents)
        # Stay constant
        rank = rank_group[name][0][2]
        taxid = rank_group[name][0][1]

        # Match rows to header
        summary.append(f"{name}\t{taxid}\t{rank}\t{percents}")

    return summary


def parse_file(fname, store={}):
    "Parse file and group its contents into a store."

    with open(fname, "r") as infile:
        for row in infile:

            # Inner loop only lasts 6 iterations and only groups the 'rank' column
            for idx, item in enumerate(row.split("\t")):
                if idx == GROUP_WITH['rank']:
                    # Append the fname to later put value in correct column
                    row = clean_row(row.split("\t"))
                    percent = [row.pop(-1) + "\t" + os.path.basename(fname)]
                    store.setdefault(item.strip(), []).append(row + percent)


def summarize_results(results):
    "Summarize result found in data_dir by grouping them."


    header = generate_header(files=results)

    store = dict()
    for item in results:
        if os.path.isfile(item):
            fname = os.path.abspath(item)
            parse_file(fname=fname, store=store)

    for x in store:
        name_store = {}
        for row in store[x]:
            # Skip the header that was parsed out of file
            if row[0] != 'name':
                name = row[0].strip()
                name_store.setdefault(name, []).append(row)

        store[x] = summarize_group(name_store, file_columns=header[3:])

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