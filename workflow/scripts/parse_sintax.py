#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
import sys


def add_lower(taxnames, lineage, ranks, i):
    """
    Add lower level labels to assignment, e.g. adds "Unclassified.GenusA"
    to an assignment that is only classified down to genus level.
    """
    missing = [r for r in ranks if not r in lineage.keys()]
    for rank in missing:
        lineage[rank] = f"Unclassified.{taxnames[-1]}"
    return lineage


def compair(lineage1, lineage2, conf1, conf2, ranks):
    """
    Compares the assignment for two read-pairs and selects either the
    most resolved one, or (if the are equally resolved) the one with a
    higher confidence value at the lowest assigned rank
    """
    l1 = [lineage1[x] for x in ranks if not lineage1[x].startswith("Unclassified.")]
    l2 = [lineage2[x] for x in ranks if not lineage2[x].startswith("Unclassified.")]
    # Compare number of matching taxlabels
    matching = set(l1).intersection(l2)
    if len(l1) > len(l2):
        return lineage1, conf1, len(matching)
    elif len(l2) > len(l1):
        return lineage2, conf2, len(matching)
    else:
        if conf1 >= conf2:
            return lineage1, conf1, len(matching)
        else:
            return lineage2, conf2, len(matching)


def replace_ranks(dataf, replacements):
    d = {}
    new_ranks = []
    for item in replacements:
        key, value = item.split("=")
        new_ranks.append(value)
        d[key] = value
    return dataf.rename(columns=d), new_ranks

def write_krona(dataf, ranks, outfile):
    if dataf.shape[0] == 0:
        sys.stderr.write("No results\n")
        with open(outfile, "w") as fhout:
            fhout.write("")
    else:
        krona_df = pd.DataFrame(dataf.loc[:, ranks].groupby(ranks).size()).reset_index()
        krona_df.loc[:, [0] + ranks].to_csv(
            outfile, index=False, sep="\t", header=False
        )


def main(args):
    rank_translator = {
        "d": "domain",
        "k": "kingdom",
        "p": "phylum",
        "c": "class",
        "o": "order",
        "f": "family",
        "g": "genus",
        "s": "species",
    }
    if not args.ranks:
        ranks = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
    else:
        ranks = args.ranks
    taxdict = {}
    confdict = {}
    cutoff = args.cutoff

    with open(args.sintax_results, "r") as fhin:
        for line in fhin:
            # For some reason sintax sometimes writes lines with only the seqid
            try:
                seqid, taxres, strand, taxlineage = line.rstrip().split("\t")
            except ValueError:
                continue
            lineage = {}
            taxnames = []
            for i, item in enumerate(taxres.split(",")):
                rank = item.split(":")[0]
                _taxname = ":".join(item.split(":")[1:])
                taxname, _conf = _taxname.split("(")
                if float(_conf.rstrip(")")) < cutoff:
                    break
                conf = float(_conf.rstrip(")"))
                taxnames.append(taxname)
                lineage[rank_translator[rank]] = taxname
            if len(taxnames) < len(ranks):
                lineage = add_lower(taxnames, lineage, ranks, i)
            if seqid in taxdict.keys():
                lineage1 = taxdict[seqid]
                conf1 = confdict[seqid]
                lineage, conf, matching = compair(lineage1, lineage, conf1, conf, ranks)
                lineage["matching"] = matching
            confdict[seqid] = float(conf)
            taxdict[seqid] = lineage
    dataf = pd.DataFrame(taxdict).T
    dataf, new_ranks = replace_ranks(dataf, args.replace_ranks)
    dataf.to_csv(sys.stdout, sep="\t")
    if args.krona:
        if len(new_ranks) > 0:
            ranks = new_ranks
        write_krona(dataf, ranks, args.krona)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("sintax_results", type=str, help="Sintax output file")
    parser.add_argument(
        "-c", "--cutoff", type=float, help="Cutoff threshold (0.8)", default=0.8
    )
    parser.add_argument(
        "-k", "--krona", type=str, help="Write krona compatible output to file"
    )
    parser.add_argument(
        "-r", "--ranks", nargs="*", help="Ranks used in the SINTAX database"
    )
    parser.add_argument("--replace_ranks", nargs="*", help="Replace rank names in final output (<rank>=<newrank>)")
    args = parser.parse_args()
    main(args)
