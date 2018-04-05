import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from random import shuffle


def join(*args):
    return os.path.abspath(os.path.join(*args))


DATA_DIR = join(os.path.dirname(__file__), "data")


def randomize_and_count(data, niter=10, subset=10):

    percent = float(subset / 100.0)

    sample_size = int(len(data) * percent)

    nunique = []

    for n in range(niter):

        shuffle(data)

        subset = data[:sample_size]

        nunique.append(len(set(subset)))

    mean_unique = 0 if not len(nunique) else sum(nunique) / len(nunique)

    return mean_unique


def generate_plot(files, niter=10, outfile=None, show=False):

    data = dict()

    # Take 10, 20, 30 ... 100% of the data.
    # note: chaining does not work with multiple files
    #subsets = chain(range(1, 30), range(20, 110, 10))
    subsets = list(range(1, 30)) + list(range(30, 110, 10))

    for fname in files:

        df = pd.DataFrame.from_csv(fname, sep='\t', header=0)
        column = df["taxID"].tolist()

        for s in subsets:
            unique_taxids = randomize_and_count(data=column, niter=niter, subset=s)
            data.setdefault(s, []).append(unique_taxids)

    legend = [os.path.basename(fname) for fname in files]
    curves = dict()

    for subset, unique_taxids in data.items():
        for idx in range(len(unique_taxids)):
            curves.setdefault(idx, []).append(unique_taxids[idx])

    for curve, label in zip(curves, legend):

        plt.plot(list(data.keys()), curves[curve], label=label)

    plt.suptitle(f"Rarefaction curve with: niter={niter}, nsamples={len(files)}")
    plt.ylabel("Number of unique species ")
    plt.xlabel("Percentage of data sampled")
    plt.ylim(ymin=0)

    outfile = join(DATA_DIR, "plot.png") or os.path.abspath(outfile)
    plt.legend()
    plt.savefig(outfile)

    if show:
        plt.show()

    return curves



def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('files', metavar='FILES', type=str, nargs='+',
                        help='a file list')

    parser.add_argument('--niter', dest='niter',
                        help="How many times to reshuffle and take subset.",
                        type=int, default=10)

    parser.add_argument('--subset', dest='subset', default=10,
                        help="Percentage of data to take out and count every iter.",
                        type=int)
    parser.add_argument('--outfile', dest='outfile', default=False, action="store_true",
                        help="Show the plot in in a GUI window.")

    parser.add_argument('--show', dest='show', default=False, action="store_true",
                        help="Show the plot in in a GUI window.")

    if len(sys.argv) == 1:
        sys.argv.extend([join(DATA_DIR, "WI-28.rep"), '--show', '--niter=100'])

    args = parser.parse_args()

    generate_plot(files=args.files, show=args.show, outfile=args.outfile, niter=args.niter)


if __name__ == '__main__':
    main()

