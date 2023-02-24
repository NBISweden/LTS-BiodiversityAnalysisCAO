#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
import os


def add_empty(empty, lines, index):
    d = {}
    for sample in empty:
        d[sample] = [0] * lines
    _counts = pd.DataFrame(d, index=index)
    return _counts


def read_counts(files):
    counts = pd.DataFrame()
    empty = []
    for f in files:
        sample = (os.path.basename(f)).replace(".sintax.parsed.tsv", "")
        df = pd.read_csv(f, sep="\t", index_col=0)
        _taxcols = list(df.columns[df.dtypes == object])
        if len(_taxcols) > 0:
            taxcols = _taxcols
        else:
            empty.append(sample)
            continue
        _counts = df.groupby(taxcols).size().reset_index()
        _counts = _counts.rename(columns={0: sample})
        index = []
        for row in _counts.iterrows():
            index.append("|".join(row[1].loc[taxcols]))
        _counts.index = index
        _counts.drop(taxcols, axis=1, inplace=True)
        counts = pd.merge(
            counts, _counts, left_index=True, right_index=True, how="outer"
        )
    namedf = pd.DataFrame([x.split("|") for x in counts.index], columns=taxcols)
    namedf.index = counts.index
    counts = pd.merge(namedf, counts, left_index=True, right_index=True)
    return counts.fillna(0), empty


def main(args):
    counts, empty = read_counts(args.files)
    if len(empty) > 0:
        _counts = add_empty(empty, counts.shape[0], counts.index)
        counts = pd.merge(counts, _counts, left_index=True, right_index=True)
    counts.to_csv(args.outfile, index=False, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("files", nargs="+", help="One or more input files")
    parser.add_argument("-o", "--outfile", type=str, help="Output file")
    args = parser.parse_args()
    main(args)
