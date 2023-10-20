#!/usr/bin/env python3

# This script outputs the number or reads for
# each taxa that align to a target genome
# and have a blast alignment to a fish, cartilagenous
# fish, mammal or bird.

# The script takes a blast_results.txt file, a filtered.sam
# file and the 'taxon_table.csv' file as input. It outputs
# two results files: one with the taxon ids and the other
# with the species ids. Note that the 'taxon_table.csv'
# file was made during the database build pipeline.

# Note the the script has been updated
# so it runs faster and so that it works with mammals
# and birds.It also requires a blast hit to be to a
# related species.


# Run as:
# python3 filter_mapping_blast_results.py blast_result.txt filtered.sam taxon_table.csv taxon_lvl_result.txt species_lvl_result.txt

#################################################
import csv
import sys
from timeit import default_timer as timer
from ete3 import NCBITaxa
import argparse
from pathlib import Path
import tempfile
import pandas as pd


def download_taxdump(tempdir):
    """
    Downloads the NCBI taxdump file to temporary directory
    """
    start = timer()
    sys.stderr.write("Downloading taxdump file\n")
    from urllib.request import urlretrieve
    import hashlib

    url = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
    md5 = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz.md5"
    filename = f"{tempdir}/taxdump.tar.gz"
    filename_md5 = f"{tempdir}/taxdump.tar.gz.md5"
    urlretrieve(url, filename)
    urlretrieve(md5, filename_md5)
    with open(filename_md5, "r") as fh:
        md5sum = fh.read().strip().split()[0]
    with open(filename, "rb") as fh:
        md5sum_file = hashlib.md5(fh.read()).hexdigest()
    if md5sum != md5sum_file:
        raise ValueError(
            f"md5sum of downloaded {filename} ({md5sum_file}) does not match {url} ({md5sum})"
        )
    end = timer()
    sys.stderr.write(f"{end - start} seconds to download taxdump file\n")
    return filename


def write_results(sp_result, tax_result, species_counts, taxid_counts):
    # Write species level counts to outfile
    sys.stderr.write(f"Writing species names counts to {sp_result}\n")
    x = csv.writer(open(sp_result, "w"))
    for key, value in species_counts.items():
        x.writerow([str(key), int(value)])

    # Write taxon level results to outfile
    sys.stderr.write(f"Writing taxon id counts to {tax_result}\n")
    w = csv.writer(open(tax_result, "w"))
    for key, value in taxid_counts.items():
        w.writerow([int(key), int(value)])


def parse_samfile(samfile):
    """
    Parse samfile and return a dictionary with the read as key and the contig as value
    """
    sys.stderr.write(f"Parsing samfile: {samfile}\n")
    start = timer()
    sam_dic = {}
    fh1 = open(samfile, "r")
    for l in fh1:
        l = l.strip().split("\t")
        r = str(l[0])
        contig = str(l[2])
        if r not in sam_dic.keys():
            sam_dic[r] = list()
            sam_dic[r].append(contig)
    fh1.close()
    end = timer()
    sys.stderr.write(f"{end - start} seconds to parse samfile\n")
    return sam_dic


def get_unique_contigs(sam_dic):
    """
    Make a list of unique contigs in the sam file. This is used to subset 'taxon_table.csv'
    """
    start = timer()
    sys.stderr.write("Getting unique contigs\n")
    contig_set = set()
    for i in sam_dic.values():
        contig_set.add(i[0])
    end = timer()
    sys.stderr.write(f"{end - start} seconds to get unique contigs\n")
    return contig_set


def group_taxon_table(df):
    return [df["species"].values[0], df["taxid"].values[0]]


def load_taxon_table(taxon_table):
    """
    Load full taxon table and groups by contig

    Returns a dataframe with contigs as index and a list ([species, taxid]) as values, e.g.:

        contig
    CAKALG010032195.1    [Helicolenus hilgendorfi, 143344]
    OMLJ01038425.1           [Gadiculus argenteus, 185737]
    OMPK01007635.1          [Myoxocephalus scorpius, 8097]
    CAKAMO010023714.1          [Sebastes mentella, 394696]
    OOFH01093767.1            [Merlangius merlangus, 8058]
    dtype: object
    """
    start = timer()
    sys.stderr.write("Loading taxon table\n")
    df = pd.read_csv(taxon_table, header=None, names=["contig", "species", "taxid"])
    dfg = df.groupby("contig").apply(group_taxon_table)
    end = timer()
    sys.stderr.write(f"{end - start} seconds to load taxon table\n")
    return dfg


def add_taxinfo(sam_dic, taxon_dic, group_dic, ncbi):
    """
    Add the species name and taxid based on the subsetted 'taxon_table.csv'
    file. Then, use the taxid to get the lineage. Add the group (e.g. mammalia,
    aves etc) and the taxid for the group.

    Desired output: {read:[contig,species,taxid,[lineage]], group, group taxid}
    """
    start = timer()
    sys.stderr.write("Adding taxonomic information\n")
    for contig, y in taxon_dic.items():
        ctg = str(contig)
        sp = str(y[0])
        txid = int(y[1])
        lineage = ncbi.get_lineage(txid)
        # Get group and group taxid (e.g. mammalia, aves etc)
        group = None
        gtxid = None
        for i in lineage:
            try:
                group = group_dic[i]
                gtxid = i
            except KeyError:
                continue
        # Update sam_dic
        for read_id, v in sam_dic.items():
            if ctg == v[0]:
                sam_dic[read_id].append(sp)
                sam_dic[read_id].append(txid)
                sam_dic[read_id].append(lineage)
                sam_dic[read_id].append(group)
                sam_dic[read_id].append(gtxid)
    end = timer()
    sys.stderr.write(f"{end - start} seconds to add taxonomic information\n")
    return sam_dic


def parse_blast_results(blast_result, taxid_col=9):
    """
    Parse blast results. Want a dictionary with read as key and a list of taxon
    ids from the blast results as values

    {read1:[taxid1,taxi2,taxid3...]}
    """
    start = timer()
    sys.stderr.write(f"Parsing blast results: {blast_result}\n")
    blast_dic = {}
    fh3 = open(blast_result, "r")
    for line in fh3:
        line = line.strip().split("\t")
        read = str(line[0])
        taxid = int(line[taxid_col])
        if read not in blast_dic.keys():
            blast_dic[read] = list()
            blast_dic[read].append(taxid)
        # Add taxid if not already present
        else:
            if taxid not in blast_dic[read]:
                blast_dic[read].append(taxid)
    fh3.close()
    end = timer()
    sys.stderr.write(f"{end - start} seconds to parse blast results\n")
    return blast_dic


def clean_blast_dict(blast_dic, sam_dic):
    """
    Remove any taxids from the blast_dic that match the taxid in the sam_dic.
    The blast hit must be to a different species. This is for cases where the
    genome assembly is also in ncbi nt.

    {read:[txid,txid,txid...]}
    """
    start = timer()
    sys.stderr.write("Cleaning blast dictionary\n")
    for x, y in blast_dic.items():
        for a, b in sam_dic.items():
            try:
                staxid = b[2]
            except IndexError:
                continue
            if staxid in y:
                y.remove(staxid)
    end = timer()
    sys.stderr.write(f"{end - start} seconds to clean blast dictionary\n")
    return blast_dic


def get_blast_group(blast_dic, group_dic, ncbi):
    """
    Go through the blast dictionary and look up the family level taxids. Go
    through the blast dictionary and make a new dictionary with the group level
    classifcation (e.g. mammalia, aves, fish etc). Note, 'gtxid' is the group
    level taxid.

    {read:[gtxid,gtxid,gtxi]}
    """
    start = timer()
    sys.stderr.write("Getting blast group\n")
    blast_group = {}
    # Get the lineage for each of the taxids
    for r, txids in blast_dic.items():
        for t in txids:
            lineage = ncbi.get_lineage(t)
            # Check if any of the values in lineage
            # are mammals, fish, birds etc
            for x in lineage:
                if x in group_dic.keys() and r not in blast_group.keys():
                    blast_group[r] = set()
                    blast_group[r].add(x)
                elif x in group_dic.keys() and r in blast_group.keys():
                    blast_group[r].add(x)
    end = timer()
    sys.stderr.write(f"{end - start} seconds to get blast group\n")
    return blast_group


def count_taxa(sam_dic, blast_group):
    """
    Check that the blast results correspond to the group (e.g. mammalia, aves
    etc) from the sam file. Go through blast_group dictionay. If value in the
    set matches the group taxid in the sam_dic, add the species to the species
    results. If, the species is already in there, add +1
    """
    start = timer()
    sys.stderr.write("Counting taxa\n")
    species_counts = {}
    taxid_counts = {}
    for r, v in sam_dic.items():
        read = r
        try:
            sp = v[1]
        except IndexError:
            continue
        txid = v[2]
        gtxid = v[5]
        # Check that the group level txids match.
        for x, y in blast_group.items():
            if x == read and gtxid in y:
                if sp not in species_counts:
                    species_counts[sp] = 1
                elif sp in species_counts:
                    species_counts[sp] += 1
                if txid not in taxid_counts:
                    taxid_counts[txid] = 1
                elif txid in taxid_counts:
                    taxid_counts[txid] += 1
    end = timer()
    sys.stderr.write(f"{end - start} seconds to count taxa\n")
    return species_counts, taxid_counts


def main(args):
    assert (
        len(args.samfile)
        == len(args.blast_result)
        == len(args.sp_result)
        == len(args.tax_result)
    ), "Number of samfiles, blast results, species results and taxid results must be the same"
    # Dictionary for birds, mammals, bony fish and
    # cartilagenous fish classification
    samfiles = sorted(args.samfile)
    blast_results = sorted(args.blast_result)
    species_results = sorted(args.sp_result)
    tax_results = sorted(args.tax_result)
    group_dic = {
        7898: "actinopterygii",
        7777: "chondrichthyes",
        40674: "mammalia",
        8782: "aves",
    }
    if not args.taxdb:
        taxdb = Path(tempfile.gettempdir()) / "taxonomy.sqlite"
        taxdb.touch()
        taxdb = str(taxdb)
        # Download taxdump file
        taxdump_file = download_taxdump(tempfile.gettempdir())
        ncbi = NCBITaxa(dbfile=taxdb, taxdump_file=taxdump_file)
    else:
        taxdb = args.taxdb
        ncbi = NCBITaxa(dbfile=taxdb)
    taxon_table = load_taxon_table(args.taxon_table)
    for i, samfile in enumerate(samfiles):
        sys.stderr.write("\n")
        blast_result = blast_results[i]
        sp_result = species_results[i]
        tax_result = tax_results[i]
        sys.stderr.write(f"Processing {samfile}\n")
        sam_dic = parse_samfile(samfile)
        contig_set = get_unique_contigs(sam_dic)
        # Get subset of taxon_table as a dictionary
        taxon_dic = taxon_table.loc[list(contig_set)].to_dict()
        sam_dic_tax = add_taxinfo(sam_dic, taxon_dic, group_dic, ncbi)
        blast_dic = parse_blast_results(blast_result, args.taxid_col)
        blast_dic = clean_blast_dict(blast_dic, sam_dic_tax)
        blast_group = get_blast_group(blast_dic, group_dic, ncbi)
        species_counts, taxid_counts = count_taxa(sam_dic_tax, blast_group)
        write_results(sp_result, tax_result, species_counts, taxid_counts)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--blast_result",
        dest="blast_result",
        nargs="+",
        help="Blast results file",
        required=True,
    )
    parser.add_argument(
        "--samfile",
        dest="samfile",
        help="Samfile",
        required=True,
        nargs="+",
    )
    parser.add_argument(
        "--taxon_table", dest="taxon_table", help="Taxon table", required=True
    )
    parser.add_argument(
        "--db",
        dest="taxdb",
        help="Taxon database",
    )
    parser.add_argument(
        "--sp_result",
        dest="sp_result",
        nargs="+",
        help="Species level results",
        required=True,
    )
    parser.add_argument(
        "--tax_result",
        dest="tax_result",
        nargs="+",
        help="Taxon id results",
        required=True,
    )
    parser.add_argument(
        "--taxid_col", dest="taxid_col", type=int, default=9, help=argparse.SUPPRESS
    )
    args = parser.parse_args()
    main(args)
