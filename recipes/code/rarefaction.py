import sys
import os
import pandas as pd
from random import shuffle

from . import plotter

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


def generate_curves(files, niter=10, outfile=None):

    data = dict()

    # Take 10, 20, 30 ... 100% of the data.
    # note: chaining does not work with multiple files
    #subsets = chain(range(1, 30), range(20, 110, 10))
    # TODO: ask about uneven buckets.
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

    title = f"Rarefaction curve with: niter={niter}, nsamples={len(files)}"
    ylabel = "Number of unique species."

    outfile = outfile if outfile.endswith(".png") else f"{outfile}.png"
    plotter.rarefactor_plot(curves=curves, legend=legend, data=data, outfile=outfile,
                            title=title, ylabel=ylabel)

    return curves



def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('files', metavar='FILES', type=str, nargs='+',
                        help='a file list')

    parser.add_argument('--niter', dest='niter',
                        help="How many times to reshuffle and take subset.",
                        type=int, default=10)
    parser.add_argument('--outfile', dest='outfile', default=False, action="store_true",
                        help="Show the plot in in a GUI window.")

    parser.add_argument('--show', dest='show', default=False, action="store_true",
                        help="Show the plot in in a GUI window.")

    if len(sys.argv) == 1:
        sys.argv.extend([join(DATA_DIR, "WI-28.rep"), '--show', '--niter=100'])

    args = parser.parse_args()

    generate_curves(files=args.files, outfile=args.outfile, niter=args.niter)



if __name__ == '__main__':
    main()

