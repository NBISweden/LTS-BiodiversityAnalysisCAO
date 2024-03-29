# LTS-BiodiversityAnalysisCAO
Code repository for the support project P_Snoeijs-Leijonmalm_2205.

## Overview

This workflow generates taxonomic and genomic profiles of metagenomic and
metatranscriptomic samples. It consists of two main tracks: _genome mapping_ and
_marker gene taxonomic assignments_. 

In the genome mapping track preprocessed reads are mapped against a set of
target species genomes followed by a filtering step using nucleotide BLAST
searches against the `nt` database. Read counts are collated and normalized to
counts per million (CPM) reads.

In the marker gene track preprocessed reads are first mapped against a database
of cytochrome oxidase subunit 1 (COI) reference sequences followed by filtering
away reads with low percent identity match or short alignment length. In
addition, a filtering step is performed to remove potential prokaryotic reads
before the remaining reads are classified taxonomically using the sintax
algorithm (implemented in `vsearch`) and a COI reference database.

## Installation

1. Clone the repository
```bash
git clone git@github.com:NBISweden/LTS-BiodiversityAnalysisCAO.git
```

2. Install the software environment
```
conda env create -f environment.yml
```

3. Activate the environment
```
conda activate biodivcao
```

## Setup & configuration

### Configuration file
The workflow behaviour can be set using an external config file in `YAML`
format. A default config file is supplied under `config/config.yml`. 

The two primary input parameters for the workflow are

- a samples file listing the samples to be analysed and
- a `mappings:` entry in the config file specifying which reference databases
  (genomes and/or marker genes) to use for mapping and taxonomic profiling

### Samples file
The samples to be used in the workflow should be specified in a comma-separated
file that you point to with the config parameter `sample_list` in the config
file, _e.g._:

```yaml
sample_list: "mysample_list.csv"
```

The samples file should have the following format:

```
sample_name,fwd_libs,rev_libs,lib_type
mysample1,path/to/mysample1_R1.fastq.gz,/path/to/mysample1_R2.fastq.gz,DNA
mysample2,path/to/mysample2_R1.fastq.gz,/path/to/mysample2_R2.fastq.gz,RNA
```

The `sample_name` column specifies sample names used in the workflow and do not
have to match the file pattern in the actual fastq files.

The `fwd_libs` and `rev_libs` columns should contain the path to fastq files for
the forward and reverse reads of each sample, respectively. If there are several
fastq files for a sample you can specify these by separating them with a
semicolon (`;`) in the `fwd_libs` and `rev_libs` columns, _e.g._:

```
sample_name,fwd_libs,rev_libs,lib_type
mysample3,path/to/mysample3-1_R1.fastq.gz;path/to/mysample3-2_R1.fastq.gz,/path/to/mysample3-1_R2.fastq.gz;path/to/mysample3-2_R1.fastq.gz,DNA
```

The `lib_type` column specifies whether the sample has reads from `DNA` or `RNA`
sequences. Currently, the workflow only maps `DNA` reads against the reference
genomes database in the [Target species track](#Target_species_track).

### Specifying reference databases

The config entry `mappings:` specifies which reference databases to use. This
config entry can have two types of nested entries: `genomes` or `marker_genes`.
The `genomes` entry specifies the reference genomes database to use in the
[Target species track](#Target_species_track) and the `marker_genes` entry
specifies the reference marker genes database to use in the [Marker gene
track](#Marker_gene_track). An example is shown below:

```yaml
mappings:
  genomes:
    target_species:
      fasta: "resources/genome_index/all_verts_whuman/target_species.fasta"
      taxon_table: "resources/genome_index/all_verts_whuman/taxon_table.csv"
      to_include: "resources/genome_index/all_verts_whuman/to_include.bed"
      mapq: 30
  marker_genes:
    coinr:
      fasta: "resources/coinr/COInr.fasta"
      ani_cutoff: 90
      min_len: 80
```

In this example, `target_species` is the name used for the genomes database. It
has three parameters: `fasta`, `taxon_table` and `to_include`. The `fasta`
parameter specifies the path to the fasta file with the reference genomes. The
`taxon_table` parameter specifies the path to the taxon table file with contigs
and their taxonomic information for genomes in the database. The `to_include`
parameter specifies the path to the bed file with the contigs to the actual
target species included in the reference genomes database. Finally, the `mapq`
parameter is the minimum map quality required for reads at the mapping stage of
the workflow. For more details see [Target species track](#target-species-track).

For the marker gene track, the example above specifies a marker gene database
named `coinr`. Each database entry has three parameters: `fasta`,
`ani_cutoff` and `min_len`. The `fasta` parameter specifies the path to the
fasta file with the reference marker genes. The `ani_cutoff` parameter specifies
the minimum average nucleotide identity (ANI) between a read and a reference
marker gene to be kept for downstream analysis. The `min_len` parameter
specifies the minimum alignment length between a read and a reference marker
gene to be kept for downstream analysis. For more details see [Marker gene
track](#marker-gene-track).

> **Note** 
Note that you can use any name for the databases. The names
`target_species` and `coinr` are just examples. You can also include any number
of database entries under the `genomes:` and `marker_genes:` entries. The
workflow will run the mapping and taxonomic profiling for each database entry.

#### Target species track
The part of the workflow that maps reads to genomes of target species (_e.g_
fish) requires a genomes database (for initial mapping) and a nucleotide BLAST
database (for filtering).  

**Target species genomes database** To identify reads originating from target
species of interest (_e.g._ fish) you must supply a fasta file with these
genomes. The relevant config parameters are:

```yaml
mappings:
  genomes:
    target_species:
      fasta: "resources/genome_index/all_verts_whuman/target_species.fasta"
      taxon_table: "resources/genome_index/all_verts_whuman/taxon_table.csv"
      to_include: "resources/genome_index/all_verts_whuman/to_include.bed"
```

The `fasta` and `taxon_table` parameters specify paths to the fasta file and
seqid->taxid mapfile, respectively. The `taxon_table` file should be a
comma-separated file with sequence ids (matching the id in the fasta file
header) in the first column, species name in the second column, species taxon id
in the third column and family taxon id in the fourth column. Examples of the two file contents are shown below:

_fasta_
```
>AESE010000001.1 Leucoraja erinacea LER_WGS_1_CONTIG_1, whole genome shotgun sequence
ACCATCATCAGCGATGATAGGGTTTACTTTAGAAATTATTCTGAGCAGAAACTTGAAGCTAAACCTCTGAAAATTTGCAA
ATCCTTTTGTATCTTTTGCAGCTGAGTTTAGAATTGTACCAAAGCACTCTAATCTTTTTTCGCAATTGATAAAACTGATG
CTTCTGTAGCACGTGGAGTAAACACTTCCTCTGAATTAAAGTCCTGTTGAAATGGTACAGAAACACTGCCAGAAGAACAC
TTGGCACCTCAATCACAATATTCTTGTAGCACTTGGTGCGTTTA
>CAKAOH010043426.1 Sebastiscus tertius genome assembly, contig: S_tertius_139529_jcf7180000782172, whole genome shotgun sequence
tcgggagcaatttggggttaagtgtcttgctcaaggacacatcgacatgtgaccggagcaaagggatcgaaccaccgacc
ttccagttgatggacgacctgctctacctacctctgagccacagtcgccCCTGAGGTACCTGAGGTAGATGAGGAtggag
gtgggggtggaggaaggGTGGAGACCTAGAGGGTGATGGGGGTGAGAGGGAGGtcgagagagaggaagggttaGTAAATg
gttattcacacacacattcacacactgatggcagaggctgccatgcaaggtgccaacctgcccat
```

_taxon\_table_
```
CM054007.1,Leucoraja erinacea,7782,30475
CAKAOH010000016.1,Sebastiscus tertius,1472224,274692
```
The `to_include` bed file that includes all of the contigs that are of interest
for downstream analysis. It could include an entry for all of the sequences in
the target species genome database or just a subset. Here, it is used to filter
alignments to human and herring. That is, the bedfile contains all of the
contigs that are included in the target species genomes database except for
human and herring.

_to\_include_
```
KB228878.1	0	17703537
KB228879.1	0	17817800
KB228880.1	0	9835568
```

_mapq_: The `mapq` parameter specifies the minimum map quality required for
reads to be kept for downstream analysis.

The target species genome database can be built in any way you like, as long as
it conforms to the file format outlined here. A GitHub repository with code to
generate such a database can be found
[here](https://github.com/NBISweden/LTS-BuildGenomesDatabaseCAO) and you're free
to use that repo as a starting point.

**BLAST database** Reads mapped to the target species database are used as
queries in a `blastn` search against a BLAST-formatted database and the results
are used as an additional filtering step to only keep reads with a high
likelihood of originating from fish genomes. Specifically, in this step reads
must have a match to sequences in the BLAST database originating from rayfinned
(actinopterygii) or cartilagenous (chondrichthyes) fish.

The relevant config parameters are:

```yaml
blast:
  db: "/sw/data/blast_databases/nt"
  evalue: 1e-30
  max_target_seqs: 200
```

The `db` parameter points to an 'alias' of the BLAST database (in this case the
actual database has several files `nt.*.nhd`, `nt.*.nhi`, `nt.*.nhr`,
`nt.*.nin`, `nt.*.nnd`, `nt.*.nni`, `nt.*.nog` and `nt.*.nsq` under
`/sw/data/blast_databases/`).

The `evalue` parameter specifies the maximum e-value of hits to keep in the
`blastn` search.

The `max_target_seqs` parameter specifies the maximum hits to keep for each
query. 

#### Marker gene track

The part of the workflow that maps reads against a marker gene database (_e.g._
the [cytochrome c oxidase subunit
I](https://en.wikipedia.org/wiki/Cytochrome_c_oxidase_subunit_I) (COI) gene)
requires:

- a fasta file of reference marker genes (for initial mapping), 
- a KrakenUniq database (for filtering) and 
- a SINTAX-formatted fasta file (for taxonomic assignments)

**Marker gene fasta file** The marker gene track of the workflow requires a
fasta file with reference sequences. The relevant config parameters are:

```yaml
mappings:
  marker_genes:
    coinr:
      fasta: "resources/coinr/COInr.fasta"
      ani_cutoff: 90
      min_len: 80
```

The `fasta` parameter specifies the path to the fasta file with the
reference marker genes. The marker gene track of the workflow starts with
mapping preprocessed reads against sequences in this fasta file. The
`ani_cutoff` and `min_len` parameters specify the minimum average nucleotide
identity (ANI) and minimum alignment length between a read and a reference
marker gene to be kept for downstream analysis. Reads passing these filters are
then passed to KrakenUniqt to filter out reads that are likely to originate from
contaminants. This second filtering can be modified by the`krakenuniq:` config
entry (see below). The remaining reads are then classified taxonomically using
the sintax algorithm (implemented in `vsearch`). 


**KrakenUniq database** 

[KrakenUniq](https://github.com/fbreitwieser/krakenuniq) is used in the workflow
primarily to filter sequences that are likely to originate from contaminants
such as prokaryotes or humans. The following config parameters control the
behaviour of this filtering:

```yaml
krakenuniq:
  db: "/sw/data/KrakenUniq_data/prebuilt/krakendb-2022-06-16-STANDARD"
  exclude: [2,2157,9604]
```

The `db` parameter specifies the path to a krakenuniq database. You can either
download prebuilt databases from
[here](https://benlangmead.github.io/aws-indexes/k2) or build a database on your
own (see the
[KrakenUniq](https://github.com/fbreitwieser/krakenuniq/blob/master/MANUAL.md))
and
[Kraken2](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)
manuals). The path specified with `db` must contain the files `database.idx`,
`database.kdb` and `database.kdb.counts`.

The `exclude` parameter is a list of taxonomic ids to exclude from downstream
analyses. These are passed to `krakenuniq-extract-reads` and only reads not
assigned to these taxonomic ids (or their children) are kept. By default, the
taxa excluded are Bacteria (taxid: 2), Archaea (taxid: 2157) and Hominidae
(taxid: 9604).  

**SINTAX database** In the marker gene mapping track of the workflow, reads are
classified taxonomically using the `sintax` algorithm implemented in
[vsearch](https://github.com/torognes/vsearch). This requires that a fasta file
is supplied and formatted in a specific way. You may use the same fasta file as
in the marker gene mapping step (described above) as long as it conforms to the
format required by `sintax` (see below).

The workflow includes steps to download a COI reference database using a
non-redundant database created using the
[mkCOInr](https://github.com/meglecz/mkCOInr/) tool. To use this database,
simply set the `fasta:` parameter under your marker gene database entry to
`resources/coinr/COInr.sintax.fasta`. Also use the same setting for `db:`
parameter under the sintax config entry (see below). The workflow will then
automatically download the database and format it for use with `sintax`.

If you want to use a custom database for the sintax classifications the fasta
file should have taxonomic information in sequence headers formatted as in the
example below:

```bash
>KJ592702_1;tax=k:Eukaryota,p:Rhodophyta,c:Florideophyceae,o:Hapalidiales,f:Mesophyllumaceae,g:Mesophyllum,s:Mesophyllum_lichenoides
```

This format is required by the sintax algorithm and allows sequences to be
classified at different ranks (here the ranks `kingdom` (k), `phylum` (p) etc.)
are used.

The behaviour of the `sintax` algorithm can be controlled using the following
config entry:

```yaml
sintax:
  cutoff: 0.8
  ranks: ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
  db: "resources/coinr/COInr.sintax.fasta"
```

The `cutoff` parameter specifies the minimum confidence cutoff used by sintax to
assign taxonomy to sequences (see the original [usearch
docs](https://drive5.com/usearch/manual/tax_conf.html)).

The `ranks` parameter specifies which ranks to use in the taxonomic
classification output. Make sure to use the same ranks as those given in the
fasta headers of the SINTAX fasta file.

## Running the workflow

When you've set up the required resources you're ready to start the workflow
with the command:

```bash
snakemake --configfile yourconfigfile.yml -j 1
```

This will start the workflow using settings supplied in `yourconfigfile.yml`
(update this to match your setup) and using at most one CPU core.

### Configuration profiles
The workflow comes with two configuration profiles: `local` and `slurm` for
running locally (_e.g._ on your laptop) or on a compute cluster with the SLURM
workload manager, respectively. You specify which of these profiles to use by
supplying either `--profile local` or `--profile slurm` to the command line
call, _e.g._:

```
snakemake --configfile yourconfigfile.yml --profile slurm
```

We strongly suggest that you use these profiles as they contain several useful
Snakemake command line settings. 

In order to use the `slurm` profile you should update the file
`slurm/config.yaml` with your SLURM account id. Modify the file by opening it in
your favourite text editor. Make sure that the `default-resources` entry has the
correct SLURM account id.

```yaml
default-resources: "slurm_account=naiss2023-1-000"
```

Now you can run the workflow as:

```bash
snakemake --profile slurm
```
