#!/usr/bin/env python

from ete3 import NCBITaxa
import sys
dbfile = sys.argv[1]
ncbi = NCBITaxa(dbfile=dbfile, taxdump_file=None)
