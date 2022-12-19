import sys, os
from os.path import join as pjoin
import shutil
from subprocess import call
from os.path import basename
import yaml
from sys import stderr
from Bio import SeqIO
import pandas
import numpy
from tempfile import mkdtemp

sys.path.append(os.getcwd())
from workflow.scripts.utils import title2log, freetxt_line

### This script is based in the mapping.py script from https://github.com/moritzbuck/0053_metasssnake2 commit 3ca68f087d7f0faa65c3400e14a9779cbb18b468

script, mapping_file, out_folder, temp_folder, threads, map_name, logfile = sys.argv

mappings_df = pandas.read_csv(mapping_file, index_col=0)
mappings = mappings_df.to_dict(orient="index")

settings = mappings[map_name]

settings.update(
    {
        "out_folder": out_folder,
        "temp_folder": temp_folder,
        "threads": threads,
    }
)

os.makedirs(temp_folder, exist_ok=True)

settings["temp_folder"] = mkdtemp(prefix=pjoin(settings["temp_folder"], "indexing_"))
settings["logfile"] = logfile
settings["envfile"] = logfile.replace("run_", "run_env")
settings["settingsfile"] = logfile.replace("run_", "run_settings")

locals().update(settings)

call(f"conda env export > {envfile}", shell=True)
with open(settingsfile, "w") as handle:
    yaml.dump(settings)


title2log("reference to temp_folder", logfile)

shutil.copy(reference, pjoin(temp_folder, "ref.fna"))

title2log(f"indexing {reference} to temp_folder", logfile)
if method == "bwa-mem2":
    call(
        f"bwa-mem2 index -p {temp_folder}/index {temp_folder}/ref.fna >> {logfile}  2>&1",
        shell=True,
    )
    with open(f"{temp_folder}/index", "w") as handle:
        pass
if method == "bowtie2":
    call(
        f"bowtie2-build --threads {threads} {temp_folder}/ref.fna {temp_folder}/index >> {logfile}  2>&1",
        shell=True,
    )
    with open(f"{temp_folder}/index", "w") as handle:
        pass
if method == "bbmap.sh":
    pass
if method == "minimap2":
    call(
        f"minimap2 -I 30G -t {threads} -d {temp_folder}/index {temp_folder}/ref.fna",
        shell=True,
    )

title2log("Cleaning up and moving the index", logfile)

for f in os.listdir(temp_folder):
    if f.startswith("index"):
        shutil.copy(
            pjoin(temp_folder, f), pjoin(out_folder, f"results/mappings/{map_name}/")
        )

shutil.rmtree(temp_folder)
