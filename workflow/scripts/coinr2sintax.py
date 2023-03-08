#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
import sys
from Bio.SeqIO import parse
from tqdm import tqdm


def read_taxa(f):
    sys.stderr.write(f"Reading taxonomy from {f}\n")
    taxdf = pd.read_csv(f, sep="\t", index_col=0, header=None,
                        names=["seq", "kingdom", "phylum", "class",
                               "order", "family", "genus", "species"],
                        dtype={"seq": str})
    df = taxdf.loc[taxdf.species != "s__"]
    sys.stderr.write(f"{df.shape[0]}/{taxdf.shape[0]} seqs with "
                     f"species assignment\n")
    return df


def main(args):
    taxdf = read_taxa(args.taxa)
    for col in tqdm(taxdf.columns, unit=" columns", desc="Reformatting columns"):
        prefix = f"{col[0]}__"
        rep = prefix[0]
        taxdf[col] = taxdf[col].str.replace(prefix,
                                            f"{rep}:").str.rstrip(";")
    sys.stderr.write("Creating lineage\n")
    lineagedf = taxdf.assign(
        lineage=taxdf["kingdom"] + "," + taxdf["phylum"] + "," +
                taxdf["class"] + "," + taxdf["order"] + "," +
                taxdf["family"] + "," + taxdf["genus"] + "," +
                taxdf["species"]

    ).loc[:, "lineage"]
    sys.stderr.write(f"Reading sequences from {args.seqs}\n")
    n = 0
    bp = 0
    with sys.stdout as fhout:
        for record in tqdm(parse(args.seqs, "fasta"), unit=" seqs"):
            try:
                header = f">{record.id};tax={lineagedf.loc[record.id]}"
                fhout.write(f"{header}\n{record.seq}\n")
                n += 1
                bp += len(record.seq)
            except KeyError:
                continue
    sys.stderr.write(f"Wrote {n} sequences, {bp} total bp\n")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("seqs", type=str, help="Sequence file")
    parser.add_argument("taxa", type=str, help="Taxonomy table")
    args = parser.parse_args()
    main(args)
