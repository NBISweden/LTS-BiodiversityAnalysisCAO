# LTS-BiodiversityAnalysisCAO
Code repository for the support project P_Snoeijs-Leijonmalm_2205.

## Overview

This workflow generates taxonomic and genomic profiles of metagenomic and 
metatranscriptomic samples. It consists of two main tracks: _genome mapping_ and
_marker gene taxonomic assignments_. 

In the genome mapping track preprocessed reads are mapped against a set of 
target species genomes followed by a filtering step using nucleotide BLAST searches
against the `nt` database. Read counts are collated and normalized counts per 
million (CPM) reads.

In the marker gene track preprocessed reads are first mapped against a database
of cytochrome oxidase subunit 1 (COI) reference sequences followed by filtering
away reads with low percent identity match or short alignment length. In addition,
a filtering step is performed to remove potential prokaryotic reads before
the remaining reads are classified taxonomically using the sintax algorithm
(implemented in `vsearch`) and the COI reference database.

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

## Testrun

```bash
snakemake --profile test -c 1
```