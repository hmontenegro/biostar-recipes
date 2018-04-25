

import csv
import os
import sys

import pandas as pd
from recipes.code import utils

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')



def map_name(name, files, delim='\t'):
    "Maps a column name to an index after checking if it is the same across files."
    idx_set = set()

    for fname in files:
        stream = csv.reader(open(fname, 'rt'), delimiter=delim)
        header = next(stream)
        idx_set.add(header.index(name))

    idx = list(idx_set)[0]
    if len(idx_set) > 1:
        msg = f"The column name does not have the same index across files."
        raise Exception(msg)

    return idx


def colnames(fnames):
    names = [os.path.basename(fname) for fname in fnames]
    names = [fname.split(".")[0] for fname in names]
    return names


def print_kreport(df, outdir=None, rank=''):

    rankmap = dict(S="Species", G="Genus", F='Family', C='Class', D='Domain')
    ranks = rank or 'SGFCD'
    pd.set_option('display.expand_frame_repr', False)

    for rank in ranks:
        subset = utils.get_subset(df, rank)
        label = rankmap.get(rank, 'Unknown')

        path = os.path.join(str(outdir), f'{label.lower()}_classification.csv')
        path = sys.stdout if not outdir else path

        subset.to_csv(index=False, path_or_buf=path)



def generate_table(allkeys, keyidx, storage, is_kreport=False):

    table = []
    for key, fields in allkeys.items():
        collect = list(reversed(fields[3:]))
        collect = collect if is_kreport else fields[0:3]

        for data in storage:
            value = data.get(key, [0]* (keyidx + 1))[keyidx]
            value = data.get(key, [0])[0] if is_kreport else value

            value = float(value)
            collect.append(value)
        table.append(collect)

    return table


def tabulate(files, rankidx=None, keyidx=4, cutoff=1, has_header=True, is_kreport=False):
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
        if has_header:
            next(stream)
        data = dict()
        for row in stream:
            data[row[keyidx]] = [elem.strip() for elem in row]
        storage.append(data)

    # Collect all taxid keys across all data.
    allkeys = {}
    for data in storage:
        allkeys.update(data)

    # The final table that can be printed and further analyzed
    table = generate_table(allkeys=allkeys, keyidx=keyidx, storage=storage,
                           is_kreport=is_kreport)

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
                        help="Name of column to combine across all files")

    parser.add_argument('--cutoff', dest='cutoff', default=0.1,
                        help="The sum of rows have to be larger than the cutoff to be registered default=%(default)s.",
                        type=float)

    parser.add_argument('--is_kreport', default=False, action='store_true',
                        help="The files to be analyzed are centrifuge kraken reports.")

    parser.add_argument('--outdir', dest='outdir', type=str,
                        help="Directory name to write data out to." )

    if len(sys.argv) == 1:
        join = lambda path: os.path.join(DATA_DIR, path)
        kreport_sample = [join('centrifuge-1.txt'), join('centrifuge-2.txt'), '--is_kreport']
        tsv_sample = [join('sample-1.tsv'), join('sample-2.tsv'), '--column=numReads', '--cutoff=0.0' ]
        sys.argv.extend(tsv_sample)

    args = parser.parse_args()

    if args.is_kreport:
        # Special case to handle kraken reports
        df = tabulate(files=args.files, cutoff=args.cutoff, rankidx=3, keyidx=4, is_kreport=True)
        print_kreport(df, outdir=args.outdir)
    else:
        # Map the column name to an id.
        assert args.column, "--column argument required."
        colid = map_name(name=args.column, files=args.files)
        df = tabulate(files=args.files, cutoff=args.cutoff, keyidx=colid)
        df.to_csv(index=False, path_or_buf=sys.stdout)




if __name__ == '__main__':
    main()
