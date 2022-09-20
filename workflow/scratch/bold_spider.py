	from tqdm  import tqdm
from os.path import join as pjoin
from Bio import SeqIO
from ete3 import NCBITaxa
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import requests
from Bio import Entrez
import xmltodict, json

Entrez.email = "murumbii@gmail.com"

query = 'ctad[All Fields] AND "bacteria"[porgn] AND (gene_nucleotide_pos[filter] AND "srcdb refseq"[Properties] AND alive[prop])'
gene_ids = Entrez.read(Entrez.esearch(db="gene", term = query , retmax=200000))['IdList']

stuff = list(Entrez.efetch(db="gene", id=list(gene_ids), retmode='xml'))
dd = xmltodict.parse("".join([b.decode().strip() for b in stuff]))
genes = dd['Entrezgene-Set']['Entrezgene']

genes = {gene['Entrezgene_gene']['Gene-ref']['Gene-ref_locus-tag'] for gene in genes}

genome_ids = Entrez.read(Entrez.esearch(db="nuccore", term = " OR ".join(genes) , retmax=1_000_000))['IdList']
genomes = [ss for ss in SeqIO.parse(Entrez.efetch(db="nuccore", id=list(genome_id), retmode = "text")]
ctdAs = [s  for genome_id in tqdm(genome_ids) for ss in SeqIO.parse(Entrez.efetch(db="nuccore", id=list(genome_id), rettype = "fasta_cds_na", retmode='text'), "fasta") if "gene=ctdA" in ss.description]
genome_dat = [Entrez.read(Entrez.esummary(db="nuccore", id= g)) for g in tqdm(genome_ids)]


def seq_complex(seq, k=2):
	twomers = re.findall(k*".", seq )
	stretches = 0
	prev_twomer = ""
	lens_twos = len(twomers)

	for i in twomers:
		if i != prev_twomer:
			stretches += 1
		prev_twomer = i


	entrops = {twomer : 0 for twomer in set(twomers) }
	for mer in twomers:
		entrops[mer] += 1/lens_twos
	entrops = sum([0 if (p == 0 or p >= 1 )else -p*log(p) - (1-p)*log(1-p) for p in entrops.values()])
	return {"stretch_complex" : stretches/lens_twos, "entropy" : entrops, "k" : k}

rank = "class"
base_url = "https://www.boldsystems.org/index.php/API_Public/sequence?taxon={taxon}"
chosen_marker = "COI-5P"
base_dir = "/home/moritz/projects/mosaic/M002_EFICA/"
taxato_dl = "Animalia"

db_path = pjoin(base_dir, "data", f"BOLD-all_{taxato_dl}.fasta")
curdb_path = pjoin(base_dir, "data", f"BOLD-clean_COI5P_{taxato_dl}.fasta")
dadadb_path = pjoin(base_dir, "data", f"BOLD-clean_COI5P_{taxato_dl}-dada.fasta")

with open("/home/moritz/dbs/bold/taxa.txt") as handle:
	taxa_todl = [l.strip() for l in handle.readlines()]

ncbi = NCBITaxa()

taxid = list(ncbi.get_name_translator(taxa_todl).values())
taxid = list({tt for t in taxid for tt in t})

all_descendents = {tt for t in taxid for tt in ncbi.get_descendant_taxa(t, intermediate_nodes = True) }
all_ranks = ncbi.get_rank(all_descendents)

filtered = [tax for tax in all_descendents if all_ranks[tax] == rank]

taxa_split = ncbi.get_taxid_translator(filtered)
taxa_split = list(ncbi.get_taxid_translator(filtered).values())

all_lines = []

for taxon in tqdm(taxa_split):
	resp = requests.get(base_url.format(taxon = taxon))
	all_lines += resp.text.split("\r\n")

with open(db_path, "w") as handle:
	handle.writelines([l + "\n" for l in all_lines])

all_seqs = [s for s in tqdm(SeqIO.parse(db_path, "fasta"))]

db_taxa = {s.description.split("|")[1] for s in all_seqs}
db_markers = {s.description.split("|")[2] for s in all_seqs}

valid_taxa = ncbi.get_name_translator(db_taxa)

valid_seqs = [s for s in all_seqs if s.description.split("|")[1] in valid_taxa and s.description.split("|")[2] == chosen_marker]
for s in valid_seqs:
	s.seq.replace("-", "")
SeqIO.write(valid_seqs, curdb_path, "fasta")

branks = [ 'kingdom',
 'phylum',
 'class',
 'order',
 'family',
 'genus',
 'species'
]
dada_refdb = []
genome2taxid = {f[0]['Caption'] : int(f[0]['TaxId']) for f in genome_dat}
for s in tqdm(ctdAs):
	dada_s = s.seq
	taxid = genome2taxid[s.id.split("|")[1].split(".")[0]]
	lineage = ncbi.get_lineage(taxid)
	taxa = ncbi.get_taxid_translator(lineage)
	ranks = ncbi.get_rank(lineage)
	taxmap = {ranks[k] : taxa[k] for k in lineage if ranks[k] in branks}
	seq_id = ";".join([taxmap.get(r, "NA").replace(" ", "_") for r in branks]).split(";NA")[0]
	dada_refdb += [SeqRecord(id = seq_id, description = "", seq = dada_s)]# s.description, seq = dada_s)]


for s in tqdm(valid_seqs):
	dada_s = s.seq
	taxid = list(ncbi.get_name_translator([s.description.split("|")[1]]).values())[0][0]
	lineage = ncbi.get_lineage(taxid)
	taxa = ncbi.get_taxid_translator(lineage)
	ranks = ncbi.get_rank(lineage)
	taxmap = {ranks[k] : taxa[k] for k in lineage if ranks[k] in branks}
	seq_id = ";".join([taxmap.get(r, "NA").replace(" ", "_") for r in branks]).split(";NA")[0]
	dada_refdb += [SeqRecord(id = seq_id, description = "", seq = dada_s)]# s.description, seq = dada_s)]

for s in ctdAs:
	dada_s = s.seq
	taxid = list(ncbi.get_name_translator([s.description.split("|")[1]]).values())[0][0]

SeqIO.write(dada_refdb, dadadb_path, "fasta")

uniqs_path = dadadb_path.replace("-dada", "-dada_uniqs")
rep_path = dadadb_path.replace("-dada", "-dada_reps")

f"usearch11 -fastx_uniques {dadadb_path} -fastaout {uniqs_path} -sizeout "
f"usearch11 -cluster_otus {uniqs_path} -otus {rep_path}"
f"sed -i 's/;size.*//' {rep_path}"

all_reps = [s for s in tqdm(SeqIO.parse(rep_path, "fasta"))]
for i,seq in enumerate(all_reps):
	seq.id = f"rep_{i}"
	seq.description = ""

SeqIO.write(all_reps, dadadb_path.replace("-dada","-mapping_ref"), "fasta")

rsnippet = """
library(dada2)
library(seqinr)

seqs = unlist(read.fasta(file = "/home/moritz/github/PhyloMagnet/test.fasta", as.string=TRUE, seqonly=TRUE))
"""
