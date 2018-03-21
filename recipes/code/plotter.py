"""
Plotter for different types of data
"""
import os
import matplotlib

# Offline mode.
HAS_DISPLAY = False

if not HAS_DISPLAY:
    # Turn off interactive display.
    matplotlib.use('Agg')

import matplotlib.pylab as plt
import numpy as np
import pandas as pd

# colidx is the column where the data starts.
# labidx is the column index for the labels
def heatmap(data, colidx=3, labidx=0, fname='heatmap.png'):

    # Based on: https://stackoverflow.com/questions/14391959/heatmap-in-matplotlib-with-pcolor

    plt.rcParams.update({'figure.autolayout': True})
    #plt.figure(figsize=(10, 5))

    df = pd.DataFrame()

    names = list(data.columns)
    label = names[labidx]
    names = [ label ] + names[colidx:]

    # A simpler dataframe with only labels and values
    for name in names:
        df[name] = data[name]

    df = df.set_index(label)

    #df = (df - df.mean()) / (df.max() - df.min())

    # Transform the scale to log.
    df = np.log(df+1)

    fig, ax = plt.subplots(figsize=(12, 6))

    heatmap = ax.pcolor(df, cmap=plt.cm.Blues, alpha=0.8)

    rnum, cnum = df.shape

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(rnum) + 0.5, minor=False)
    ax.set_xticks(np.arange(cnum) + 0.5, minor=False)

    #ax.invert_yaxis()
    ax.xaxis.tick_top()

    labels = list(df.columns)
    ax.set_xticklabels(labels, minor=False)
    ax.set_yticklabels(df.index, minor=False)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    plt.xticks(rotation=90)

    plt.savefig(f'{fname}')

    if HAS_DISPLAY:
        # Pop a window in non-offline mode.
        plt.show()



# colidx is the column where the data starts.
def horizontal_bars(data, colidx=3, labidx=0, fname='plot.png'):
    # Set plotting parameters.
    plt.rcParams.update({'figure.autolayout': True})
    plt.figure(figsize=(20, 10), dpi=100)

    # How rows and columns in the dataframe.
    rnum, cnum = data.shape

    # The tick positions.
    ypos = np.arange(rnum)

    # The total width for bars.
    space = 0.3

    # How many bars will there be.
    n = cnum - 3

    # Width of one bar.
    width = (1 - space) / n

    # Shift the bar to center
    shift = (n - 1) * width / 2

    # Generate a barplot for each column past the third.
    for i in range(n):
        label = data.columns[colidx + i]
        values = data[label]
        npos = ypos + i * width - shift
        plt.barh(npos, values, width, label=label)

    # Finish up the plot.
    plt.legend()
    # The label column.
    labels = data[data.columns[labidx]]
    plt.yticks(range(len(labels)), labels)
    plt.title(f'Read Classification')
    plt.xlabel("Percent reads")
    plt.savefig(f'{fname}')

    if HAS_DISPLAY:
        # Pop a window in non-offline mode.
        plt.show()


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('input', metavar='FILES', type=str, nargs=1,
                        help='The file to be plotted')

    parser.add_argument('--type', dest='type', default='hbar',
                        help="Plot type",
                        type=str)

    parser.add_argument('--output', dest='output', default='plot.png',
                        help="Plot file name",
                        type=str)


if __name__ == '__main__':
    main()