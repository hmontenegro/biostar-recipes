

import csv
import os
import sys
import warnings

import pandas as pd
from recipes.code import utils

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')



def is_digit(string):
    "Check if a string holds an int or float."

    try:
        float(string)
    except ValueError:
        return False
    else:
        return True


def map_name(name, files, delim='\t'):
    "Maps a column name to an index after checking if it is the same across files."
    idx_set = set()

    for fname in files:
        stream = csv.reader(open(fname, 'rt'), delimiter=delim)
        header = next(stream)
        idx_set.add(header.index(name))

    idx = list(idx_set)[0]
    if len(idx_set) > 1:
        msg = f"The column name does not have the same index across files. Using:{idx}"
        warnings.warn(msg, SyntaxWarning)

    return idx


def colnames(fnames):
    names = [os.path.basename(fname) for fname in fnames]
    names = [fname.split(".")[0] for fname in names]
    return names


def print_csv(df, outdir, combined_col=None):
    if os.path.exists(outdir):
        name = os.path.join(outdir, f'combined_{combined_col}.csv')
        df.to_csv(path_or_buf=name, index=False)
    else:
        print(df.to_csv(index=False))


def print_data(df, rank='', outdir=None, by_rank=False, combined_col=None):

    if not by_rank:
        print_csv(df=df, outdir=outdir, combined_col=combined_col)
        return

    rankmap = dict(S="Species", G="Genus", F='Family', C='Class', D='Domain')
    ranks = rank or 'SGFCD'
    pd.set_option('display.expand_frame_repr', False)

    for rank in ranks:
        subset = utils.get_subset(df, rank)
        label = rankmap.get(rank, 'Unknown')

        if os.path.exists(outdir):
            path = os.path.join(outdir, f'{label.lower()}_classification.csv')
            subset.to_csv(index=False, path_or_buf=path)
        else:
            print(f'#Rank: {label}')
            print(subset.to_csv(index=False))


def tabulate(files, rankidx=None, keyidx=4, cutoff=1, select=None):
    """
    Summarize result found in data_dir by grouping them.
    """

    # Collect all data into a dictionary keyed by keyID
    storage = []
    for fname in files:
        stream = csv.reader(open(fname, 'rt'), delimiter="\t")
        # Keep only known rank codes
        if rankidx != None:
            stream = filter(lambda x: x[rankidx] != '-', stream)
        data = dict()
        for row in stream:
            data[row[keyidx]] = [elem.strip() for elem in row]
        storage.append(data)

    # Collect all taxid keys across all data.
    allkeys = {}
    for data in storage:
        allkeys.update(data)

    # The final table that can be printed and further analyzed
    table = []
    for key, fields in allkeys.items():
        if not is_digit(string=key):
            continue

        collect = list(reversed(fields[3:])) if select != None else fields[0:3]
        for data in storage:
            value = data.get(key, [0]* (keyidx + 1))[keyidx]
            value = value if select == None else data.get(key, [0])[select]

            value = float(value)
            collect.append(value)
        table.append(collect)

    # Filter table by cutoffs
    cond1 = lambda row: sum(row[3:]) > cutoff
    table = list(filter(cond1, table))

    # Sort by reverse of the abundance.
    compare = lambda row: (row[2], sum(row[3:]))
    table = sorted(table, key=compare, reverse=True)

    # Make a panda dataframe
    columns = ["name", "taxid", "rank"] + colnames(files)
    df = pd.DataFrame(table, columns=columns)

    # Attempt to fill in common names at species level.
    fname = "/export/refs/alias/fishalias.txt"
    df = utils.alias(df=df, fname=fname, key1='name', key2='sciname', name='name')

    return df


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('files', metavar='FILES', type=str, nargs='+',
                        help='a file list')

    parser.add_argument('--column', dest='column', type=str,
                        help="Name of column to combine across all files",
                        )

    parser.add_argument('--cutoff', dest='cutoff', default=0.1,
                        help="The sum of rows have to be larger than the cutoff to be registered default=%(default)s.",
                        type=float)

    parser.add_argument('--outdir', dest='outdir', default='combine_report',
                        help="Directory name to write data out to.",
                        type=str)
    if len(sys.argv) == 1:
        full = lambda path: os.path.join(DATA_DIR, path)
        txt_sample = [full('centrifuge-1.txt'), full('centrifuge-2.txt')]
        tsv_sample = [full('1000-MiniXFish.tsv'), full('WI-10.tsv'),
                      '--column=numReads', '--cutoff=0.0', ]

        sys.argv.extend(tsv_sample)

    args = parser.parse_args()

    if not args.column:
        df = tabulate(files=args.files, cutoff=args.cutoff, rankidx=3, select=0)
        # Prints by rank, and produces one file for each rank
        by_rank = True
    else:
        # Map the column name to an id.
        colid = map_name(name=args.column, files=args.files)
        df = tabulate(files=args.files, cutoff=args.cutoff, keyidx=colid)
        # Print one combined file
        by_rank = False

    print_data(df, outdir=args.outdir, combined_col=args.column, by_rank=by_rank)

if __name__ == '__main__':
    main()
