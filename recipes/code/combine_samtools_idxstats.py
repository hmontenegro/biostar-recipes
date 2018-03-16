"""

This program is used to process idx stats files created by samtools idxstats.

"""

import os
import sys

import pandas as pd

pd.set_option('display.expand_frame_repr', False)
pd.options.display.float_format = '{:.1f}'.format


def get_colnames(fnames):
    names = [os.path.basename(fname) for fname in fnames]
    names = [fname.split(".")[0] for fname in names]
    return names


def plot(data, args):
    from recipes.code import plotter

    # Plot a heatmap
    plotter.heatmap(data=data, colidx=2, fname="heatmap.png")


# Prints a dataframe
def print_data(df):
    print(df.to_csv(index=False, float_format='%.1f'))


ACCESSION, SIZE, MAPPED, UNMAPPED, SUM = ['accession', 'size', 'mapped', 'unmapped', 'sum']


def tabulate(args):
    """
    Summarize the index stats.
    """

    files = sorted(args.files)
    colnames = [ACCESSION, SIZE, MAPPED, UNMAPPED]
    res = None
    for fname in files:
        df = pd.read_table(fname, header=None)
        df.columns = colnames
        rnum, cnum = df.shape
        unmapped = df[UNMAPPED].sum()
        df.loc[rnum-1, MAPPED] = unmapped
        values = df[MAPPED]
        percent = values / values.sum() * 100
        df[MAPPED] = percent
        if res is None:
            res = df[[ACCESSION, SIZE, MAPPED]]
        else:
            df = df[[ACCESSION, MAPPED]]
            res = pd.merge(res, df[[ACCESSION, MAPPED, ]], on=ACCESSION, how='inner')

    # Rename the columns
    fnames = get_colnames(files)

    colnames = [ACCESSION, SIZE] + fnames
    res.columns = colnames

    # Sort table by the sum of columns
    res[SUM] = res[fnames].sum(axis=1)
    res = res.sort_values(by=SUM, ascending=False)
    res = res[res[SUM] > args.cutoff]

    # Drop the sum column
    res = res.drop(SUM, axis=1)

    return res


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('files', metavar='FILES', type=str, nargs='+',
                        help='a file list')

    parser.add_argument('--cutoff', dest='cutoff', default=0.1,
                        help="The sum of rows have to be larger than the cutoff to be registered default=%(default)s.",
                        type=float)

    if len(sys.argv) == 1:
        sys.argv.extend(['data/1000-MiFish_R1.fq.idxstats.txt', 'data/1000-MiFish_R2.fq.idxstats.txt'])

    args = parser.parse_args()

    df = tabulate(args)

    # Print the data to screen.
    print_data(df)

    plot(df, args=args)


if __name__ == '__main__':
    main()
