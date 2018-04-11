import os
<<<<<<< HEAD
import pandas as pd
from random import shuffle

from . import plotter
=======
import sys
from random import shuffle

import matplotlib.pyplot as plt
import pandas as pd

>>>>>>> 961a6850920510340ba1129af9de36095c762464

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


def generate_curves(files, niter=10, outfile=None):

    # Percents of the data to compute the counts over
    percents = list(range(0, 20, 5)) + list(range(20, 120, 20))

<<<<<<< HEAD
    # Take 10, 20, 30 ... 100% of the data.
    # note: chaining does not work with multiple files
    #subsets = chain(range(1, 30), range(20, 110, 10))
    # TODO: ask about uneven buckets.
    subsets = list(range(1, 30)) + list(range(30, 110, 10))
=======
    plt.figure(figsize=(12, 8))
>>>>>>> 961a6850920510340ba1129af9de36095c762464

    for fname in files:
        df = pd.read_csv(fname, sep='\t', header=0)
        column = df["taxID"].tolist()

        x = percents
        y = [randomize_and_count(data=column, niter=niter, percent=p) for p in percents]
        label = os.path.basename(fname)

<<<<<<< HEAD
    title = f"Rarefaction curve with: niter={niter}, nsamples={len(files)}"
    ylabel = "Number of unique species."

    outfile = outfile if outfile.endswith(".png") else f"{outfile}.png"
    plotter.rarefactor_plot(curves=curves, legend=legend, data=data, outfile=outfile,
                            title=title, ylabel=ylabel)
=======
        plt.plot(x, y, 'bo-', label=label)


    plt.suptitle(f"Rarefaction curve with: niter={niter}, nsamples={len(files)}")
    plt.ylabel("Number of unique species ")
    plt.xlabel("Percentage of data sampled")
    plt.ylim(ymin=0)

    plt.legend()

    plt.savefig(outfile)

    if show:
        plt.show()
>>>>>>> 961a6850920510340ba1129af9de36095c762464

    return (x, y)


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('files', metavar='FILES', type=str, nargs='+',
                        help='a file list')

    parser.add_argument('--niter', dest='niter',
                        help="How many times to reshuffle and take subset.",
<<<<<<< HEAD
                        type=int, default=10)
    parser.add_argument('--outfile', dest='outfile', default="plot.png", type=str,
                        help="Show the plot in in a GUI window.")
=======
                        type=int, default=150)

    parser.add_argument('--plot', dest='outfile', default="rarefaction.png",
                        help="The name of the plot file.")
>>>>>>> 961a6850920510340ba1129af9de36095c762464

    parser.add_argument('--show', dest='show', default=False, action="store_true",
                        help="Show the plot in in a GUI window.")

    if len(sys.argv) == 1:
        sys.argv.extend([join(DATA_DIR, "WI-28.rep"), f'--outfile={join(DATA_DIR, "test.png")}', '--niter=100'])

    args = parser.parse_args()

    generate_curves(files=args.files, outfile=args.outfile, niter=args.niter)



if __name__ == '__main__':
    main()
