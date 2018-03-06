import os





def parse_file(rep_file, store={}):
    "Parse a rep file and populate dict with contents "

    # What do the columns mean? to help when grouping


    print(rep_file)
    return



def findfiles(location, collect=[], exts=(".txt",)):
    "Walk a directory and only return files with desired extension. "


    for item in os.scandir(location):

        if item.is_dir():
            findfiles(item.path, collect=collect)
        else:
            # Check if file is in list of extensions.
            has_ext = True in list(map(lambda ext: item.path.endswith(ext), exts))
            if has_ext and exts:
                collect.append(os.path.abspath(item.path))


    return collect




def summarize_results(data_dir, group_by):
    "Summarize result found in data_dir. Exclusively looks at .rep and .rep.txt files."

    store = dict()

    all_files = findfiles(location=data_dir, exts=(".rep", ".rep.txt"))

    for fname in all_files:

        parse_file(rep_file=fname, store=store)


    # Order


    return




def main():
    import sys
    from argparse import ArgumentParser

    parse = ArgumentParser()
    parse.add_argument('--dir', dest='indir', required=True,
                       help="""Directory containing .rep/.rep.txt files to be parsed.""",
                       type=str)
    parse.add_argument('--group_by', dest='group_by', default="domain",
                       help="""Group resulting report in specific manner.""",
                       # Domain is always on the fourth column and species on the last?
                       type=str, choices={"domain": 3, "species": 5})

    args = parse.parse_args()

    if len(sys.argv) == 1:
        sys.argv.append('-h')

    if os.path.isdir(args.indir):
        summarize_results(data_dir=args.indir, group_by=args.group_by)
    else:
        parse.print_help()
        print('--dir needs to be a valid directory.')
        exit()



if __name__ == '__main__':

    import sys

    try:
        main()
    except (BrokenPipeError, KeyboardInterrupt) as e:
        sys.stderr.close()