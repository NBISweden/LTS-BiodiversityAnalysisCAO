# LTS-BiodiversityAnalysisCAO
Code repository for the support project P_Snoeijs-Leijonmalm_2205.

## Overview

This workflow generates taxonomic and genomic profiles of metagenomic and 
metatranscriptomic samples. It consists of two main tracks: _genome mapping_ and
_marker gene taxonomic assignments_. 

In the genome mapping track preprocessed reads are mapped against a set of 
target species genomes followed by a filtering step using nucleotide BLAST searches
against the `nt` database. Read counts are collated and normalized to counts per 
million (CPM) reads.

In the marker gene track preprocessed reads are first mapped against a database
of cytochrome oxidase subunit 1 (COI) reference sequences followed by filtering
away reads with low percent identity match or short alignment length. In addition,
a filtering step is performed to remove potential prokaryotic reads before
the remaining reads are classified taxonomically using the sintax algorithm
(implemented in `vsearch`) and a COI reference database.

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

### Samples file
The samples to be used in the workflow should be specified in a 
comma-separated file that you point to with the config parameter 
`sample_list` in the config file, _e.g._:

```yaml
sample_list: "mysample_list.csv"
```

The samples file should have the following format:

```
sample_name,fwd_libs,rev_libs,lib_type
mysample1,path/to/mysample1_R1.fastq.gz,/path/to/mysample1_R2.fastq.gz,DNA
mysample2,path/to/mysample2_R1.fastq.gz,/path/to/mysample2_R2.fastq.gz,RNA
```

The `sample_name` column specifies sample names used in the workflow and do 
not have to match the file pattern in the actual fastq files.

The `fwd_libs` and `rev_libs` columns should contain the path to fastq files 
for the forward and reverse reads of each sample, respectively. If there are 
several fastq files for a sample you can specify these by separating them 
with a semicolon (`;`) in the `fwd_libs` and `rev_libs` columns, _e.g._:

```
sample_name,fwd_libs,rev_libs,lib_type
mysample3,path/to/mysample3-1_R1.fastq.gz;path/to/mysample3-2_R1.fastq.gz,/path/to/mysample3-1_R2.fastq.gz;path/to/mysample3-2_R1.fastq.gz,DNA
```

The `lib_type` column specifies whether the sample has reads from `DNA` or 
`RNA` sequences. Currently, the workflow only maps `DNA` reads against the 
reference genomes database in the [Target species track](#Target_species_track).

### Installing resources

The workflow depends on some external resources that should be set up prior 
to running analyses.

#### Target species track
The part of the workflow that maps reads to genomes of target species (_e.g_ 
fish) requires a genomes database (for initial mapping) and a nucleotide 
BLAST database (for filtering).  

**Target species genomes database**
To identify reads originating from target species of interest (_e.g._ fish) you 
must supply a fasta file with these genomes. The relevant config parameters 
are:

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
header) in the first column, species name in the second column and taxonomy 
ids in the third column. Examples of the two file contents are shown below:

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

_taxon_table_
```
AESE010000001.1,Leucoraja erinacea,7782
CAKAOH010043426.1,Sebastiscus tertius,1472224
```
The `to_include` bed file that includes all of the contigs that are of interest 
for downstream analysis. It could include an entry for all of the sequences in
the target species genome database or just a subset. Here, it is used to filter
alignments to human and herring. That is, the bedfile contains all of the contigs
that are included in the target species genomes database except for human and
herring.

_to_include_
```
KB228878.1	0	17703537
KB228879.1	0	17817800
KB228880.1	0	9835568
```

The target species genome database can be built in any way you like, as long 
as it conforms to the file format outlined here. A GitHub repository with 
code to generate such a database can be found 
[here](https://github.com/NBISweden/LTS-BuildGenomesDatabaseCAO) and you're free 
to use that repo as a starting point.

**BLAST database**
Reads mapped to the target species database are used as queries in a 
`blastn` search against a BLAST-formatted database and the results are used 
as an additional filtering step to only keep reads with a high likelihood of 
originating from fish genomes. Specifically, in this step reads must have a 
match to sequences in the BLAST database originating from rayfinned
(actinopterygii) or cartilagenous (chondrichthyes) fish.

The relevant config parameters are:

```yaml
blast:
  db: "/sw/data/blast_databases/nt"
  evalue: 1e-30
  max_target_seqs: 200
```

The `db` parameter points to an 'alias' of the BLAST database (in this case 
the actual database has several files `nt.*.nhd`, `nt.*.nhi`, 
`nt.*.nhr`, `nt.*.nin`, `nt.*.nnd`, `nt.*.nni`, `nt.*.nog` and `nt.*.nsq` 
under `/sw/data/blast_databases/`).

The `evalue` parameter specifies the maximum e-value of hits to keep in the 
`blastn` search.

The `max_target_seqs` parameter specifies the maximum hits to keep for each 
query. 

#### Marker gene track

The part of the workflow that maps reads against a marker gene database (_e.g._ 
the [cytochrome c oxidase subunit I](https://en.wikipedia.org/wiki/Cytochrome_c_oxidase_subunit_I)
(COI) gene) requires:

- a fasta file of reference marker genes (for initial mapping), 
- a KrakenUniq database (for filtering) and 
- a SINTAX-formatted fasta file (for taxonomic assignments)

**Marker gene fasta file**
The initial mapping of reads is done against a fasta file specified in a 
separate comma-separated file with the format:

```
map_name,reference,ani_cutoff,min_len
coinr,resources/coinr/COInr.fasta,90,80
```

You point the workflow to this file with the config parameter 
`mappings_list:`, _e.g._:

```yaml
mappings_list: "config/mappings.csv"
```

The `map_name` column in the file gives the database a name and can be set 
to anything you like. The `reference` column points to the actual fasta 
file on your system. The `ani_cutoff` and `min_len` parameters specify the 
minimum average nucleotide identity and minimum alignment length that reads in 
your samples must have to at least one of the reference sequences in the fasta 
file in order to be kept for downstream analysis.

You can add several lines to the `mappings_list` file to map reads against 
more than one fasta file. Just make sure to use different values in the 
`map_name` column for each line. You may use the same path in the 
`reference` column and simply update the `map_name`, `ani_cutoff` and 
`min_len` fields to run the initial mapping with different parameters.

**KrakenUniq database**
[KrakenUniq](https://github.com/fbreitwieser/krakenuniq) is used in the 
workflow primarily to filter sequences that are likely to originate from 
contaminants such as prokaryotes or humans. The following config parameters 
control the behaviour of this filtering:

```yaml
krakenuniq:
  db: "/sw/data/KrakenUniq_data/prebuilt/krakendb-2022-06-16-STANDARD"
  exclude: [2,2157,9604]
```

The `db` parameter specifies the path to a krakenuniq database. You can 
either download prebuilt databases from 
[here](https://benlangmead.github.io/aws-indexes/k2) or build a database on 
your own (see the [KrakenUniq](https://github.com/fbreitwieser/krakenuniq/blob/master/MANUAL.md))
and [Kraken2](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)
manuals). The path specified with `db` must contain the files `database.idx`,
`database.kdb` and `database.kdb.counts`.

The `exclude` parameter is a list of taxonomic ids to exclude from 
downstream analyses. These are passed to `krakenuniq-extract-reads` and only 
reads not assigned to these taxonomic ids (or their children) are kept. By 
default, the taxa excluded are Bacteria (taxid: 2), Archaea (taxid: 2157) and 
Hominidae (taxid: 9604).  

**SINTAX database**
In the marker gene mapping track of the workflow, reads are classified 
taxonomically using the `sintax` algorithm implemented in 
[vsearch](https://github.com/torognes/vsearch). This requires that a 
sintax-formatted fasta file is available for use. The relevant config 
parameters are:

```yaml
sintax:
  cutoff: 0.8
  ranks: ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
  replace_ranks: ["kingdom=kingdom", "phylum=phylum", "class=class", "order=order", "family=family", "genus=genus", "species=species"]
  db: "resources/coinr/COInr.sintax.fasta"
```

Here the `db` parameter is a path to a fast file with headers formatted as:

```bash
>KJ592702_1;tax=k:Eukaryota_2759_kingdom,p:Rhodophyta_2763,c:Florideophyceae_2806,o:Hapalidiales_1705611,f:Mesophyllumaceae_2784734,g:Mesophyllum_48973,s:Mesophyllum_lichenoides_1000564
```

This format is required by the sintax algorithm and allows sequences to be 
classified at different ranks (here the ranks `kingdom` (k), `phylum` (p) 
etc.) are used.

The `cutoff` parameter specifies the minimum confidence cutoff used by sintax 
to assign taxonomy to sequences (see the original 
[usearch docs](https://drive5.com/usearch/manual/tax_conf.html)).

The `ranks` parameter specifies which ranks to use in the taxonomic 
classification output. We suggest to use the same ranks as those given in 
the fasta headers of the SINTAX fasta file.

The `replace_ranks` parameter allows you to remap ranks from those specified 
in the SINTAX fasta headers. For example if you have:

```
>PLYAO072-20;tax=d:Animalia,k:Arthropoda,p:Insecta,c:Lepidoptera,o:Gelechiidae,f:Stegasta,g:Stegasta bosqueella,s:BOLD:AAA1027
```

in your fasta file but want to rename the shorthand `s` value to `bold_id` 
in the parsed sintax output you can use:

```yaml
sintax:
  replace_ranks: ["domain=kingdom", "kingdom=phylum", "phylum=class", "class=order", "order=family", "family=genus", "genus=species", "species=bold_id"]
```

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
running locally (_e.g._ on your laptop) or on a compute cluster with the 
SLURM workload manager, respectively. You specify which of these profiles to 
use by supplying either `--profile local` or `--profile slurm` to the 
command line call, _e.g._:

```
snakemake --configfile yourconfigfile.yml --profile slurm
```

We strongly suggest that you use these profiles as they contain 
several useful Snakemake command line settings. 

In order to use the `slurm` profile you should update the file 
`slurm/settings.json` with your SLURM account id. Modify the file by opening 
it in your favourite text editor:
```json
{
    "SBATCH_DEFAULTS": "account=SLURM_ACCOUNT no-requeue exclusive",
    "CLUSTER_NAME": "",
}
```

Change `SLURM_ACCOUNT` to match your SLURM account id.

When you run with the `slurm` profile the `-j` flag in the snakemake command 
specifies number of jobs that can be run in parallel. If you specify 
`--profile slurm` you can omit the `-j` flag since it is set automatically 
by the profile (default: 500).
