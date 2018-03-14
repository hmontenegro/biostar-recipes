"""

This program is used to process .tsv report files that are outputted by centrifuge-kreport.
It groups the results according to rank code and shows the abundance for each file.

The input files are expected to have this header:

    name	taxID	taxRank	genomeSize	numReads	numUniqueReads	abundance

Usage:
    
    $   python centrifuge.py --files=data/*.tsv

    name	                    taxID	taxRank	    SRR1972972_1.tsv	SRR1972972_2.tsv    ... 
    Azorhizobium caulinodans	7	    species	    0.0	                0.0                 
    Cellulomonas gilvus	        11	    species	    -                   0.0
    Pelobacter carbinolicus	    19	    species	    0.0                 -
    ...
                   
"""

import os
import sys
import glob
from numpy import arange
import matplotlib.pyplot as plt



HEADER = dict(name=0, taxID=1, taxRank=2, genomeSize=3,
              numReads=4,numUniqueReads=5,abundance=6)


def generate_header(files):

    header = ["name", "taxID", 'taxRank'] + [os.path.basename(p) for p in files]

    return header


def compute_sum(values):

    try:
         return sum(float(x) for x in filter(lambda x: x not in ('-', '0', '0.0'), values))
    except ValueError:
        return 1


def plot(data, outfile, header, title='Title', width=0.35, opacity=0.8):


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



def summarize_group(rank_group, file_columns, col_idx=-1, drop_zeros=False):
    "Summarize a nested list 'rank_group' by flattening."

    summary = []
    col_to_idx = {x:y for y,x in enumerate(file_columns)}

    for name in rank_group:
        # Start off with empty values and substitute
        # every index that has a value
        values = ["-" for _ in range(len(file_columns))]
        current = [x[col_idx] for x in rank_group[name]]
        for p in current:
            # Variable 'col' is the column name value belongs to
            value, col = p.split('\t')[0], p.split('\t')[1]
            # Substitute the '-' in correct index with 'percent'
            values[col_to_idx[col]] = value
        is_zero = not compute_sum(values=values)
        if drop_zeros and is_zero:
            continue
        rank = rank_group[name][0][2].split('\t')[0]
        taxid = rank_group[name][0][1].split('\t')[0]
        name, values = name.split('\t')[0], '\t'.join(values)

        # Match rows to header
        summary.append(f"{name}\t{taxid}\t{rank}\t{values}")

    return summary


def parse_file(fname, store={}, col_idx=HEADER['abundance']):
    "Parse file and group its contents into a store."

    with open(fname, "r") as infile:
        for row in infile:

            # Inner loop only lasts 6 iterations and only groups the 'rank' column
            for idx, item in enumerate(row.split("\t")):
                if idx == HEADER['taxRank']:
                    row = list(map(lambda val: val.strip(),row.split("\t")))
                    # Add the column name to later put value in correct place
                    row[col_idx] = row[col_idx] + "\t" + os.path.basename(fname)
                    store.setdefault(item.strip(), []).append(row)


def summarize_results(results, filter_by='', drop_zeros=False, summarize='abundance'):
    "Summarize result found in data_dir by grouping them."

    header = generate_header(files=results)
    summ_column = HEADER[summarize]

    store = dict()
    for item in results:
        if os.path.isfile(item):
            fname = os.path.abspath(item)
            # Parse each file and populate store with grouped ranks
            parse_file(fname=fname, store=store, col_idx=summ_column)

    for rank in store:
        # Group rows by name within each rank
        # then use grouped dictionary to summarize abundance for each name.
        name_store = {}
        for row in store[rank]:
            # Skip the header that was parsed out of incoming files
            if row[0].strip() != header[0] and rank not in header:
                name = row[0].strip()
                name_store.setdefault(name, []).append(row)
        # Flatten data for each name
        store[rank] = summarize_group(name_store, file_columns=header[3:],
                                      col_idx=summ_column, drop_zeros=drop_zeros)
    if filter_by:
        return dict(filter_by=store.get(filter_by, []))

    return store


def main():
    from argparse import ArgumentParser

    parse = ArgumentParser()
    parse.add_argument('--files', dest='files', required=True,
                       help="""Directory containing files to be parsed.""",
                       type=str)
    parse.add_argument('--outfile', dest='outfile',
                       help="""Output file to write summary report to.""",
                       type=str)
    parse.add_argument('--drop_zeros', dest='drop_zeros', default=False,
                       help="""Drop rows that have numerical values equal zero across all samples.""",
                       action='store_true')
    parse.add_argument('--filter', dest='filter',
                       help="""Filter summary report to one rank.""",
                       type=str)
    parse.add_argument('--plot', dest='plot', default=False,
                        help="""Produce a bar plot for each rank in summary report""",
                        action='store_true')
    parse.add_argument('--prefix', dest='prefix', default='plot',
                       help="""Prefix to the plot name.""",
                       type=str)
    parse.add_argument('--summarize', dest='summarize', default='abundance',
                       help="""Column to combine across all samples and put in summary report.""",
                       type=str,choices=("genomeSize", "numReads" ,"numUniqueReads","abundance"))

    args = parse.parse_args()
    if len(sys.argv) == 1:
        sys.argv.append('-h')

    files = glob.glob(args.files)
    summary = summarize_results(results=files, filter_by=args.filter,
                                drop_zeros=args.drop_zeros, summarize=args.summarize)
    if not args.outfile:
        print('\t'.join(generate_header(files=files)))
        for rank, values in summary.items():
            if values:
                #print('\n'.join(values))
                pass
    else:
        with open(args.outfile, "w") as outfile:
            outfile.write('\t'.join(generate_header(files=files)))
            for rank, values in summary.items():
                if values:
                    outfile.write('\n'.join(values) + "\n")
    if args.plot:
        for rank, values in summary.items():
            # Generate a plot for each rank
            plot(data=values, title=f'{rank}',outfile=args.prefix + '.png',
                 header=generate_header(files))


if __name__ == '__main__':

    try:
        main()
    except (BrokenPipeError, KeyboardInterrupt) as e:
        sys.stderr.close()