#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
from tqdm import tqdm
import sys


def concat_cols(df, ranks):
    d = df.loc[df.name, ranks].to_dict()
    return ",".join([f"{k[0]}:{d[k]}" for k in d.keys()])


def main(args):
    tqdm.pandas(desc="Making sintax headers", unit=" seq")
    dtype={
        "seqID": str,
        "taxon": str, 
        "taxID": int, 
        "taxlevel": float, 
        "superkingdom": str, 
        "superkingdom_taxID": str,
        "kingdom": str,
        "kingdom_taxID": str,
        "phylum": str,
        "phylum_taxID": str,
        "class": str,
        "class_taxID": str,
        "order": str,
        "order_taxID": str,
        "family": str,
        "family_taxID": str,
        "genus": str,
        "genus_taxID": str,
        "species": str,
        "species_taxID": str,
        "sequence": str,}
    sys.stderr.write(f"Reading input file {args.infile}\n")
    df = pd.read_csv(args.infile, sep="\t", header=0, index_col=0, dtype=dtype)
    sys.stderr.write(f"Keeping only {args.ranks} ranks\n")
    dataframe = df.loc[:, args.ranks+["sequence"]]
    # remove sequences with no taxonomic assignment at species level
    sys.stderr.write(f"Removing sequences with no species assignment\n")
    dataframe = dataframe.loc[dataframe.species==dataframe.species]
    sys.stderr.write(f"{dataframe.shape[0]} sequences remaining\n")
    strings = dataframe.groupby(level=0).progress_apply(concat_cols, ranks=args.ranks)
    seq_df = pd.merge(pd.DataFrame(dataframe.loc[:, "sequence"]), pd.DataFrame(strings, columns=["string"]), left_index=True, right_index=True)
    with open(args.outfile, 'w') as fhout:
        for row in tqdm(seq_df.iterrows(), desc="Writing output", unit=" seq"):
            fhout.write(f">{row[1].name};tax={row[1].string}\n{row[1].sequence}\n")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("infile", type=str, help="Full tsv file from format_db.pl")
    parser.add_argument("outfile", type=str, help="Sintax formatted fasta file")
    parser.add_argument("--ranks", nargs="+", default=["kingdom", "phylum", "class", "order", "family", "genus", "species"], help="Taxonomic ranks to include in the header")
    args = parser.parse_args()
    main(args)
