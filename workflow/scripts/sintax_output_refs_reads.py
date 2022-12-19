#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd


def write_reads(outfile, l):
    with open(outfile, "w") as fhout:
        for r in l:
            fhout.write(f"^{r}\n")

def main(args):
    taxdf = pd.read_csv(args.tax, sep="\t", index_col=0)
    # Get rank names
    sp_hits = taxdf.loc[taxdf[args.rank].str.startswith("Unclassified")]
    queries = list(sp_hits.index)
    sp = sp_hits[args.rank].unique()
    dbinfo = pd.read_csv(args.dbinfo, sep="\t", index_col=0)
    # Get database records
    records = dbinfo.loc[dbinfo[args.rank].isin(sp)].index
    write_reads(records, args.refs_out)
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
