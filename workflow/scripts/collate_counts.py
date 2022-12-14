#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
import os


def read_files(input, strip):
    df = pd.DataFrame()
    for f in input:
        sample = os.path.basename(f).rstrip(strip)
        try:
            _df = pd.read_csv(f, header=None, index_col=0)
            _df.columns = [sample]
        except pd.errors.EmptyDataError:
            _df = pd.DataFrame(columns=[sample])
        df = pd.merge(df, _df, left_index=True, right_index=True, how="outer")
    df.index.name = "tax"
    return df.fillna(0)


def main(args):
    counts = read_files(args.input, args.strip)
    counts.to_csv(args.outfile, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("input", nargs="+", help="Input file with counts")
    parser.add_argument("outfile", type=str, help="Write reults to outfile")
    parser.add_argument(
        "--strip", type=str, help="Strip this text from basename to get sample name"
    )
    args = parser.parse_args()
    main(args)
