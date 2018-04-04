"""

This program is used to process .tsv report files that creatred via centrifuge-kreport.

It groups the results according to rank code and produces the abundance for each file.

"""

import csv
import os
import sys

import pandas as pd
from recipes.code import utils

def colnames(fnames):
    names = [os.path.basename(fname) for fname in fnames]
    names = [fname.split(".")[0] for fname in names]
    return names


def plot(df, args):
    from recipes.code import plotter

    # Subst to species.
    data = get_subset(df, 'S')

    # Plot a heatmap
    plotter.heatmap(data=data, fname="heatmap.png")

def get_subset(df, rank=''):
    indices = df['rank'] == rank
    subset = df[indices] if rank else df
    return subset


# Prints a dataframe at a rank.
def print_data(df, rank=''):
    rankmap = dict(S="Species", G="Genus", F='Family', C='Class', D='Domain')
    ranks = rank or 'SGFCD'
    pd.set_option('display.expand_frame_repr', False)
    for rank in ranks:
        subset = get_subset(df, rank)
        label = rankmap.get(rank, 'Unknown')
        print(f'###\n### Rank: {label}\n###')
        print(subset.to_csv(index=False))


def tabulate(files, rank='', rankidx=3, keyidx=4, cutoff=1):
    "Summarize result found in data_dir by grouping them."

    # Collect all data into a dictionary keyed by keyID
    storage = []
    for fname in files:
        stream = csv.reader(open(fname, 'rt'), delimiter="\t")

        # Keep only known rank codes
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
        collect = list(reversed(fields[3:]))
        for data in storage:
            value = data.get(key, [0])[0]
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

    # Attempt to fill in common names at species level:
    fname = "/export/refs/alias/fishalias.txt"
    df = utils.alias(df=df, fname=fname, left='name', right='sciname', column='name')

    return df


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('files', metavar='FILES', type=str, nargs='+',
                        help='a file list')

    parser.add_argument('--rank', dest='rank',
                        help="Filter summary report to a taxonomic rank default.",
                        type=str, default='')

    parser.add_argument('--cutoff', dest='cutoff', default=0.1,
                        help="The sum of rows have to be larger than the cutoff to be registered default=%(default)s.",
                        type=float)

    parser.add_argument('--show', dest='show', default=False, action="store_true",
                        help="Show the plot in in a GUI window.")

    if len(sys.argv) == 1:
        sys.argv.extend(['--show', 'data/centrifuge-1.txt', 'data/centrifuge-2.txt'])

    args = parser.parse_args()

    df = tabulate(files=args.files, rank=args.rank,
                  cutoff=args.cutoff)

    # Print the data to screen.
    print_data(df)

    plot(df=df, args=args)


if __name__ == '__main__':
    main()