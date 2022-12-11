# This script outputs the number or reads for
# each taxa that have a map to a target genome
# and have a blast alignment to a fish/cartilagenous
# fish. It takes the filtered samfile, a blast result,
# and the taxon information for the contigs in the
# mapping database as input. The output is a file
# with the taxid (or species id) and a count for each

# The script take a blast_results.txt file, a filtered.sam
# file and the 'taxon_table.csv' file as input. It outputs
# two results files: one with the taxon ids and the other
# with the species ids. Note that the 'taxon_table.csv'
# file was made during the database build pipeline.

# Run as:
# python3 filter_aln_and_blast_results.py blast_result.txt filtered.sam taxon_table.csv taxon_lvl_result.txt species_lvl_result.txt

#################################################

import csv
import sys
from ete3 import NCBITaxa

ncbi = NCBITaxa()

# Infiles/outfiles
blast_result = sys.argv[1]
samfile = sys.argv[2]
taxon_table = sys.argv[3]
tax_result = sys.argv[4]
sp_result = sys.argv[5]


# Taxids for rayfinned fish and sharks/rays group
actinopterygii = int(7898)
chondrichthyes = int(7777)

# Parse blast results. Want a dictionary with
# read as key and a list of taxon ids as values
# {read1:[taxid1,taxi2,taxid3]}

blast_dic = {}

fh = open(blast_result, "r")
for line in fh:
    line = line.strip().split("\t")
    read = str(line[0])
    taxid = int(line[11])
    if read not in blast_dic.keys():
        blast_dic[read] = list()
        blast_dic[read].append(taxid)
    # Add taxid if not already present
    else:
        if taxid not in blast_dic[read]:
            blast_dic[read].append(taxid)
fh.close()

# Go through the dictionary, get the full
# lineage for each taxon id, then check if
# it matches an actinopterygii species or
# a chondrichthyes species

fish_reads = list()
shark_reads = list()
for k, v in blast_dic.items():
    list_values = v
    for i in list_values:
        lineage = ncbi.get_lineage(i)
        for x in lineage:
            if x == actinopterygii and k not in fish_reads:
                fish_reads.append(k)
            elif x == chondrichthyes and k not in shark_reads:
                shark_reads.append(k)


# Next is to filter the sam file and select
# the rows that are either 'fish' or sharks
# {read:[contig]}
results_dic = {}
fh2 = open(samfile, "r")
for l in fh2:
    l = l.strip().split("\t")
    r = str(l[0])
    contig = str(l[2])
    if r in fish_reads and r not in results_dic.keys():
        results_dic[r] = list()
        results_dic[r].append(contig)
    elif r in shark_reads and r not in results_dic.keys():
        results_dic[r] = list()
        results_dic[r].append(contig)
fh2.close()


# Next is to add the taxid and species name
# based on the contig id. Use the taxid file
# from the database build pipeline. The
# file is called 'taxon_table.csv'
# Desired output:
# {read:[contig,taxid,species]}

fh3 = open(taxon_table, "r")
for row in fh3:
    row = row.strip().split(",")
    ctg = str(row[0])
    sp = str(row[1])
    txid = int(row[2])
    for k, v in results_dic.items():
        if ctg == v[0]:
            results_dic[k].append(txid)
            results_dic[k].append(sp)
fh3.close()

# This prints a count for each taxon id
taxid_counts = {}
for i in results_dic.values():
    tid = int(i[1])
    if tid not in taxid_counts:
        taxid_counts[tid] = 1
    else:
        taxid_counts[tid] += 1

# write results to outfile
w = csv.writer(open(tax_result, "w"))
for key, value in taxid_counts.items():
    w.writerow([int(key), int(value)])

# This prints a count for each species
species_counts = {}
for j in results_dic.values():
    sid = str(j[2])
    if sid not in species_counts:
        species_counts[sid] = 1
    else:
        species_counts[sid] += 1

# Write species level counts to outfile
x = csv.writer(open(sp_result, "w"))
for key, value in species_counts.items():
    x.writerow([str(key), int(value)])
