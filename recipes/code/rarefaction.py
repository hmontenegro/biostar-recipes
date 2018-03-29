import sys
import pandas as pd
from random import shuffle





def randomize_and_count(data, niter=10, subset=10):

    percent = float(subset / 100.0)

    sample_size = len(data) * percent

    nunique = []

    for n in range(niter):

        shuffle(data)

        subset = data[:int(sample_size)]

        nunique.append(len(set(subset)))

    mean_unique = 0 if not len(nunique) else sum(nunique) / len(nunique)

    return mean_unique




def rarefactor_plot(data, show=False, outfile=None):

    print(data)


    return





def generate_curves(files, niter=10, subset=10, show=True):


    curves = dict()
    subsets = range(10, 110, subset)

    for fname in files:

        df = pd.DataFrame.from_csv(fname, sep='\t', header=0)

        for s in subsets:

            unique_taxids = randomize_and_count(data=df["taxID"].tolist(), niter=niter, subset=s)

            curves.setdefault(s, []).append(unique_taxids)

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
        sys.argv.extend(['--h'])

    args = parser.parse_args()

    curves = generate_curves(files=args.files, niter=args.niter, subset=args.subset)

    rarefactor_plot(data=curves, show=args.show, outfile=args.outfile)


if __name__ == '__main__':
    main()









