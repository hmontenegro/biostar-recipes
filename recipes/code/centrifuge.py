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

import sys, os, csv

def plot(data, outfile, header, title='Title', width=0.35, opacity=0.8):

    import matplotlib.pyplot as plt

    objects = {}
    for row in data:
        name = row.split('\t')[0].strip()
        row = row.split('\t')[3:]
        objects.setdefault(name, []).append(row)

    y_pos = arange(len(objects))
    ax, fig = plt.subplot(), plt.gcf()

    for idx, x_val in enumerate(objects.items()):
        #TODO: needs to match y_pos shape
        d = [float(x) for y in x_val[1] for x in y]
        print(d, y_pos, len(d), len(y_pos), x_val)
        # Shift to accommodate more than one bar
        # max number of bars = len (colors) = 7
        #if color > 0:
        y_pos = y_pos + width
        plt.barh(y_pos + width, d, width,
                alpha=opacity,
                 # Label and color need to come from filename.
                #color=colors[color],
                #label=n if not EASYMAP.get(n) else EASYMAP[n],
                align = 'center')

    plt.yticks(y_pos + width, objects)
    #plt.xlabel(f'{x_label}', weight='bold')
    plt.title(f'{title}', fontsize=10, weight='bold')
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles, labels, loc='right', bbox_to_anchor=(0.15, -0.1),
                  prop={'size': 8})
    plt.tight_layout()
    #fig.set_size_inches(9, 5.5)
    plt.show()
    #plt.savefig(outfile)
    print(f"Plot saved to :{os.path.abspath(outfile)}")
    print(objects)
    print(objects)
    1/0

    return


def tabulate(files, rank, column, keyid='taxID', cutoff=1):
    "Summarize result found in data_dir by grouping them."

    # Find the headers
    fname = files[0]
    stream = csv.reader(open(fname, 'rt'), delimiter="\t")
    header = next(stream)

    # Find the rank index
    rankidx = header.index(rank) if rank in header else -1

    # Collect all data into a dictionary keyed by keyID
    storage = []
    for fname in files:
        stream = csv.DictReader(open(fname, 'rt'), delimiter="\t")
        data = dict()
        for row in stream:
            data[row[keyid]] = row
        storage.append(data)

    # Collect all keys across all data.
    allkeys = {}
    for data in storage:
        allkeys.update(data)

    # The final table that can be printed and further analyzed
    table = []
    for key, fields in allkeys.items():
        collect = [ fields[name] for name in header[:4]]
        for data in storage:
            fields = data.get(key, {})
            value = int(fields.get(column, 0))
            collect.append(value)
        table.append(collect)

    # Filter table by cutoffs
    cond1 = lambda row: sum(row[4:]) > cutoff
    table = list(filter(cond1, table))

    # Filter by ranks if required
    cond2 = lambda row: (row[rankidx]) == rank
    table = list(filter(cond2, table)) if rankidx > -1 else table

    # Print the table
    basenames = [os.path.basename(fname) for fname in files]
    basenames = [os.path.splitext(fname)[0] for fname in basenames]
    outheader = header[:4] + basenames
    print("\t".join(outheader))
    for row in table:
        row = map(str, row)
        print ("\t".join(row))

    return table


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('files', metavar='FILES', type=str, nargs='+',
                        help='a file list')

    parser.add_argument('--rank', dest='rank',
                       help="Filter summary report to a taxonomic rank default.",
                       type=str, default='')

    parser.add_argument('--column', dest='column', default='numReads',
                       help="Column to combine across all samples and put in summary report default=%(default)s.",
                       type=str)

    parser.add_argument('--cutoff', dest='cutoff', default=10,
                        help="The sum of rows have to be larger than the cutoff to be registered default=%(default)s.",
                        type=int)

    parser.add_argument('--plot', dest='plot', default=False,
                        help="Produce a bar plot for each rank in summary report.",
                        action='store_true')

    parser.add_argument('--prefix', dest='prefix', default='plot',
                        help="Prefix to the plot name.",
                        type=str)

    if len(sys.argv) == 1:
        sys.argv.extend( [ 'data/SRR1972972_1.tsv', 'data/SRR1972972_2.tsv' ])

    args = parser.parse_args()


    results = tabulate(files=args.files, rank=args.rank,
                       column=args.column, cutoff=args.cutoff)

    '''
    if args.plot:
        for rank, values in summary.items():
            # Generate a plot for each rank
            plot(data=values, title=f'{rank}',outfile=args.prefix + '.png',
                 header=generate_header(files))
    '''

if __name__ == '__main__':
    main()
