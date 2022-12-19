#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
import sys


def write_reads(l, outfile, suffix=""):
    with open(outfile, "w") as fhout:
        for r in l:
            fhout.write(f"^{r}{suffix}\n")


def main(args):
    taxdf = pd.read_csv(args.tax, sep="\t", index_col=0)
    sys.stderr.write(f"Read {taxdf.shape[0]} result lines from {args.tax}\n")
    # Get rank names
    sp_hits = taxdf.loc[~taxdf[args.rank].str.startswith("Unclassified")]
    sys.stderr.write(f"{len(sp_hits)} reads assigned at {args.rank}\n")
    queries = list(sp_hits.index)
    sp = sp_hits[args.rank].unique()
    sys.stderr.write(f"{len(sp)} unique {args.rank}\n")
    sys.stderr.write(f"Reading info from {args.dbinfo}\n")
    dbinfo = pd.read_csv(args.dbinfo, sep="\t", index_col=0)
    # Get database records
    records = dbinfo.loc[dbinfo[args.rank].isin(sp)].index
    sys.stderr.write(f"{len(records)} reference records matched\n")
    write_reads(records, args.refs_out, suffix=";")
    write_reads(queries, args.queries_out)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("tax", type=str, help="Sintax parsed results")
    parser.add_argument("dbinfo", type=str, help="COIDB info file")
    parser.add_argument("refs_out", type=str, help="Outfile with reference ids")
    parser.add_argument("queries_out", type=str, help="Outfile with query ids")
    parser.add_argument(
        "--rank", type=str, default="species", help="Rank at which to filter"
    )
    args = parser.parse_args()
    main(args)
