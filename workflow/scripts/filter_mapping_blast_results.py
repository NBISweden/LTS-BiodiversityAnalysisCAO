
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
#python3 filter_mapping_blast_results.py blast_result.txt filtered.sam taxon_table.csv taxon_lvl_result.txt species_lvl_result.txt

#################################################
import csv
import sys
from ete3 import NCBITaxa

# Infiles/outfiles
blast_result = sys.argv[1]
samfile = sys.argv[2]
taxon_table = sys.argv[3]
tax_result = sys.argv[4]
sp_result = sys.argv[5]

# For running within the main workflow, use:
taxdb = sys.argv[6]
ncbi = NCBITaxa(dbfile=taxdb)
ncbi.update_taxonomy_database()

# For testing outside the main workflow, use:
#ncbi = NCBITaxa()

# Dictionary for birds, mammals, bony fish and
# cartilagenous fish classification
group_dic = {7898:'actinopterygii',7777:'chondrichthyes', 40674:'mammalia', 8782:'aves'}


# Parse samfile
# {read:[contig]}

sam_dic={}
fh1=open(samfile,"r")
for l in fh1:
	l = l.strip().split('\t')
	r = str(l[0])
	contig = str(l[2])
	if r not in sam_dic.keys():
		sam_dic[r] = list()
		sam_dic[r].append(contig)
fh1.close()


# Make a list of unique contigs in the sam file. This is
# used to subset 'taxon_table.csv'
contig_set = set()
for i in sam_dic.values():
	contig_set.add(i[0])


# Load taxon table into a dictionary based on the contigs
# in the unique list. {contig:[species,taxid]}
taxon_dic = {}
fh2=open(taxon_table,'r')
for row in fh2:
	row = row.strip().split(',')
	ctg = str(row[0])
	sp = str(row[1])
	txid = int(row[2])
	l = [sp,txid]
	if ctg in contig_set:
		taxon_dic[ctg] = l

fh2.close()

# Add the species name and taxid based on the 
# subsetted 'taxon_table.csv' file. Then, use
# the taxid to get the lineage. Add the group
# (e.g. mammalia, aves etc) and the taxid for
# the group.   
# Desired output: {read:[contig,species,taxid,[lineage]], group, group taxid}

for x,y in taxon_dic.items():
	ctg = str(x)
	sp = str(y[0])
	txid = int(y[1])
	lineage = ncbi.get_lineage(txid)
	# Get group and group taxid (e.g. mammalia, aves etc)
	for i in lineage:
		for t in group_dic.keys():
			if i == t:
				gtxid = t
				group = group_dic[t]
	# Update sam_dic
	for k,v in sam_dic.items():
		if ctg == v[0]:
			sam_dic[k].append(sp)
			sam_dic[k].append(txid)
			sam_dic[k].append(lineage)
			sam_dic[k].append(group)
			sam_dic[k].append(gtxid)

		

# Parse blast results. Want a dictionary with
# read as key and a list of taxon ids from the
# blast results as values
# {read1:[taxid1,taxi2,taxid3...]}

blast_dic = {}

fh3=open(blast_result,"r")
for line in fh3:
	line = line.strip().split('\t')
	read = str(line[0])
	taxid = int(line[9])
	if read not in blast_dic.keys():
		blast_dic[read] = list()
		blast_dic[read].append(taxid)
	# Add taxid if not already present
	else:
		if taxid not in blast_dic[read]:
			blast_dic[read].append(taxid)
fh3.close()


# Remove any taxids from the blast_dic that
# match the taxid in the sam_dic. The blast
# hit must be to a different species. This is
# for cases where the genome assembly is also
# in ncbi nt.
# {read:[txid,txid,txid...]}

for x,y in blast_dic.items():
	for a,b in sam_dic.items():
		staxid = b[2]
		if staxid in y:
			y.remove(staxid)


# Go through the blast dictionary and look up
# the family level taxids.
# Go through the blast dictionary and make
# a new dictionary with the group level classifcation
# (e.g. mammalia, aves, fish etc). Note, 'gtxid' is 
# the group level taxid {read:[gtxid,gtxid,gtxi]}

blast_group = {}

# Get the lineage for each of the taxids
for r,txids in blast_dic.items():
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



# Check that the blast results correspond to the
# group (e.g. mammalia, aves etc) from the sam file.
# Go through blast_group dictionay. If value in the set
# matches the group taxid in the sam_dic, add
# the species to the species results. If, the
# species is already in there, add +1

species_counts = {}
taxid_counts = {}

for r,v in sam_dic.items():
	read = r
	sp = v[1]
	txid = v[2]
	gtxid = v[5]
	# Check that the group level txids match.
	for x,y in blast_group.items():
		if x == read and gtxid in y:
			if sp not in species_counts:
				species_counts[sp] = 1
			elif sp in species_counts:
				species_counts[sp] += 1	
			if txid not in taxid_counts:
				taxid_counts[txid] = 1
			elif txid in taxid_counts:
				taxid_counts[txid] += 1


# Write species level counts to outfile
x = csv.writer(open(sp_result, "w"))
for key,value in species_counts.items():
	x.writerow([str(key),int(value)])


# Write taxon level results to outfile
w = csv.writer(open(tax_result, "w"))
for key,value in taxid_counts.items():
	w.writerow([int(key),int(value)])

