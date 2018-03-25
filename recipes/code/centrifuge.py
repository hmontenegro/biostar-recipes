
"""
This program is used to process .tsv report files that are outputted by centrifuge-kreport.
It groups the results according to rank code and shows the abundance for each file.
The input files are expected to have this header:
    name	taxID	taxRank	genomeSize	numReads	numUniqueReads	abundance
Usage:
    
    $   python centrifuge.py data/*.tsv
    name	                    taxID	taxRank	    SRR1972972_1.tsv	SRR1972972_2.tsv    ... 
    Azorhizobium caulinodans	7	    species	    0.0	                0.0                 
    Cellulomonas gilvus	        11	    species	    -                   0.0
    Pelobacter carbinolicus	    19	    species	    0.0                 -
    ...
                   
"""

import csv
import os
import sys

import pandas as pd


def colnames(fnames):
    names = [os.path.basename(fname) for fname in fnames]
    names = [fname.split(".")[0] for fname in names]
    return names


def plot(df, args):
    import matplotlib.pylab as plt
    import numpy as np


    plt.title(f'FOO', fontsize=10, weight='bold')
    data = get_subset(df, 'S')
    values = data['1000-MiFish_R1']
    labels = data['name']
    ypos = np.arange(len(values))


    fig, ax = plt.subplots()
    ax.barh(ypos, values)
    ax.set_yticklabels(labels)


    plt.yticks(np.arange(len(labels)))



    plt.show()


    return


def get_subset(df, rank=''):
    indices = df['rank'] == rank
    subset = df[indices] if rank else df
    return subset

# Prints a dataframe at a rank.
def print_data(df, rank=''):
    ranks = rank or 'SGFCD'
    for rank in ranks:
        subset = get_subset(df, rank)
        print(subset)
        print('-' * 80)

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
    table = sorted(table, key=compare, reverse=False)

    # Make a panda dataframe
    columns = ["name", "taxid", "rank"] + colnames(files)
    df = pd.DataFrame(table, columns=columns)

    return df


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('files', metavar='FILES', type=str, nargs='+',
                        help='a file list')

    parser.add_argument('--rank', dest='rank',
                        help="Filter summary report to a taxonomic rank default.",
                        type=str, default='')

    parser.add_argument('--cutoff', dest='cutoff', default=1,
                        help="The sum of rows have to be larger than the cutoff to be registered default=%(default)s.",
                        type=int)

    parser.add_argument('--plot', dest='plot', default=False,
                        help="Produce a bar plot for each rank in summary report.",
                        action='store_true')

    parser.add_argument('--prefix', dest='prefix', default='plot',
                        help="Prefix to the plot name.",
                        type=str)

    if len(sys.argv) == 1:
        sys.argv.extend(['--plot', 'data/1000-MiFish_R1.fq.txt', 'data/1000-MiFish_R2.fq.txt'])

    args = parser.parse_args()

    df = tabulate(files=args.files, rank=args.rank,
                     cutoff=args.cutoff)

    # Print the data to screen.
    print_data(df)

    if args.plot:
        plot(df=df, args=args)


if __name__ == '__main__':
    main()