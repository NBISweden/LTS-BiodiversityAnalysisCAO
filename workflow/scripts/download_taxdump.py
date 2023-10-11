#!/usr/bin/env python

from ete3 import NCBITaxa
dbfile = "resources/taxonomy/taxonomy.sqlite"
taxdump_file = "resources/taxonomy/taxdump.tar.gz"
ncbi = NCBITaxa(dbfile=dbfile, taxdump_file=None)
