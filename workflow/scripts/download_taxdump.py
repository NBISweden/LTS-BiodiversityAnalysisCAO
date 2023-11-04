#!/usr/bin/env python

from ete3 import NCBITaxa
import sys
from pathlib import Path

dbfile = sys.argv[1]
Path(dbfile).touch()

ncbi = NCBITaxa(dbfile=dbfile, taxdump_file=None)
