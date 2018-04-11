import os
import sys
from random import shuffle

import matplotlib.pyplot as plt
import pandas as pd


def join(*args):
    return os.path.abspath(os.path.join(*args))


DATA_DIR = join(os.path.dirname(__file__), "data")


def randomize_and_count(data, niter=10, percent=10):
    # Turn the percent into fraction.
    frac = float(percent / 100.0)

    # The selection size.
    sample_size = int(len(data) * frac)

    counts = []

    for n in range(niter):
        shuffle(data)
        subset = data[:sample_size]
        counts.append(len(set(subset)))

    mean_count = 0 if not len(counts) else sum(counts) / len(counts)

    return mean_count


def generate_plot(files, niter=10, outfile=None, show=False):

    # Percents of the data to compute the counts over
    percents = list(range(0, 20, 5)) + list(range(20, 120, 20))

    plt.figure(figsize=(12, 8))

    for fname in files:
        df = pd.read_csv(fname, sep='\t', header=0)
        column = df["taxID"].tolist()

        x = percents
        y = [randomize_and_count(data=column, niter=niter, percent=p) for p in percents]
        label = os.path.basename(fname)

        plt.plot(x, y, 'bo-', label=label)


    plt.suptitle(f"Rarefaction curve with: niter={niter}, nsamples={len(files)}")
    plt.ylabel("Number of unique species ")
    plt.xlabel("Percentage of data sampled")
    plt.ylim(ymin=0)

    plt.legend()

    plt.savefig(outfile)

    if show:
        plt.show()

    return (x, y)


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('files', metavar='FILES', type=str, nargs='+',
                        help='a file list')

    parser.add_argument('--niter', dest='niter',
                        help="How many times to reshuffle and take subset.",
                        type=int, default=150)

    parser.add_argument('--plot', dest='outfile', default="rarefaction.png",
                        help="The name of the plot file.")

    parser.add_argument('--show', dest='show', default=False, action="store_true",
                        help="Show the plot in in a GUI window.")

    if len(sys.argv) == 1:
        sys.argv.extend([join(DATA_DIR, "WI-28.rep"), '--show', '--niter=100'])

    args = parser.parse_args()

    generate_plot(files=args.files, show=args.show, outfile=args.outfile, niter=args.niter)


if __name__ == '__main__':
    main()