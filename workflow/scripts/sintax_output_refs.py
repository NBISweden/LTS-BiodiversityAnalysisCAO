#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd


def main(args):
    taxdf = pd.read_csv(args.tax, sep="\t", index_col=0)
    # Get rank names
    sp = taxdf.loc[taxdf[args.rank].str.startswith("Unclassified")][args.rank].unique()
    dbinfo = pd.read_csv(args.dbinfo, sep="\t", index_col=0)
    # Get database records
    records = dbinfo.loc[dbinfo[args.rank].isin(sp)].index
    with open(args.outfile, "w") as fhout:
        for r in records:
            fhout.write(f"{r}\n")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("tax", type=str, help="Sintax parsed results")
    parser.add_argument("dbinfo", type=str, help="COIDB info file")
    parser.add_argument("outfile", type=str, help="Outfile with reads")
    parser.add_argument(
        "--rank", type=str, default="species", help="Rank at which to filter"
    )
    args = parser.parse_args()
    main(args)
