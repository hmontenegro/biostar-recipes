import pandas as pd


def alias(df, fname, left='name', right='name', column=''):
    """
    Replaces the 'lect' column of a dataframe with the `column` value where 'right' matches `left`.
    """

    try:

        frame = pd.read_table(fname, delimiter="\t")

        alias = dict()
        for name1, name2 in zip(frame[right], frame[column]):
            alias[name1] = name2

        names = []
        for value in df[left]:
            names.append(alias.get(value, value))
        df[left] = names

    except Exception as exc:
        print(f"Aliasing not performed: {exc}")
        pass

    return df
