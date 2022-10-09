import sys, os
from os.path import join as pjoin
import shutil
from subprocess import call
from os.path import basename
import yaml
from datetime import datetime
from math import floor, ceil
from sys import stderr
from Bio import SeqIO
import pandas
import numpy
from tempfile import mkdtemp
sys.path.append(os.getcwd())

from workflow.scripts.utils import title2log, freetxt_line

### This script is based in the mapping.py script from https://github.com/moritzbuck/0053_metasssnake2 commit 3ca68f087d7f0faa65c3400e14a9779cbb18b468


script , sample_file, mapping_file, out_folder, temp_folder, threads, sample, map_name, logfile = sys.argv

sample_df = pandas.read_csv(sample_file, index_col=0)
samples = sample_df.to_dict(orient="index")

mappings_df = pandas.read_csv(mapping_file, index_col=0)
mappings = mappings_df.to_dict(orient="index")


settings = mappings[map_name]
settings.update(samples[sample])

settings.update({
    "out_folder" :  out_folder,
    "temp_folder" : temp_folder,
    "threads" : threads,
})


settings['temp_folder'] = mkdtemp(prefix = pjoin(settings['temp_folder'], "coi_mapping_"))
os.makedirs(temp_folder, exist_ok=True)

settings['logfile'] = logfile
settings['envfile'] = logfile.replace("run_", "run_env")
settings['settingsfile'] = logfile.replace("run_", "run_settings")

settings['fwd_libs'] = settings['fwd_libs'].split(";")
settings['rev_libs'] = settings['rev_libs'].split(";")

os.makedirs(os.path.dirname(settings['logfile']), exist_ok=True)

locals().update(settings)

call(f"conda env export > {envfile}", shell=True)
with open(settingsfile, "w") as handle:
    yaml.dump(settings)


os.makedirs(temp_folder, exist_ok=True)

title2log("copying index to temp_folder", logfile)

for f in os.listdir( pjoin(out_folder, f"results/mappings/{map_name}/")):
    if f.startswith("index"):
        shutil.copy(pjoin(out_folder, f"results/mappings/{map_name}/", f), pjoin(temp_folder))

freetxt_line(f"Starting to copy {len(fwd_libs)} pairs of fastqs", logfile)

for ff, rr in zip(fwd_libs,rev_libs):
    title2log(f"copying {ff} and {rr} to temp_folder", logfile)
    call(f"""
    unpigz -kc {ff}> {temp_folder}/fwd.fastq 2>> {logfile}
    unpigz -kc {rr} > {temp_folder}/rev.fastq 2>> {logfile}
    """, shell=True)

call(f"fastp -y -l {min_len} -h /dev/null -j {out_folder}/fastp.json  --in1 {temp_folder}/fwd.fastq --in2 {temp_folder}/rev.fastq --out1 {temp_folder}/clean_fwd.fastq --out2 {temp_folder}/clean_rev.fastq -w {threads}  >> {logfile} 2>&1", shell= True)


freetxt_line("Starting mappings", logfile)
title2log(f"mapping {sample} to ref with {method}", logfile)
if method == "bwa-mem2":
    call(f"""
    bwa-mem2 mem -t {threads} {temp_folder}/index  -o {temp_folder}/mapping.sam {temp_folder}/clean_fwd.fastq {temp_folder}/clean_rev.fastq 2>> {logfile}
    samtools view  -b -S -@{threads}  {temp_folder}/mapping.sam | samtools sort -@ 24 -o {temp_folder}/mapping_pairs.bam - >> {logfile} 2>&1
    rm {temp_folder}/mapping.sam 2>> {logfile}
    """, shell=True)
if method == "bowtie2" :
    call(f"""
    bowtie2 -p {threads} -x  {temp_folder}/index --very-sensitive -S {temp_folder}/mapping.sam  -1 {temp_folder}/clean_fwd.fastq -2 {temp_folder}/clean_rev.fastq >> {logfile} 2>&1
    samtools view  -b -S -@{threads}  {temp_folder}/mapping.sam | samtools sort -@ 24 -o {temp_folder}/mapping_pairs.bam - >> {logfile} 2>&1
    bowtie2 -p {threads} -x  {temp_folder}/index --very-sensitive -S {temp_folder}/mapping.sam  -U {temp_folder}/unp.fastq >> {logfile} 2>&1
    samtools view -b -S -@{threads}  {temp_folder}/mapping.sam | samtools sort -@ 24 -o {temp_folder}/mapping_unpaired.bam - >> {logfile} 2>&1
    rm {temp_folder}/mapping.sam 2>> {logfile}
    samtools merge -f -t {threads} {temp_folder}/mapping.bam  {temp_folder}/mapping_pairs.bam  {temp_folder}/mapping_unpaired.bam >> {logfile}  2>&1
    rm {temp_folder}/mapping_pairs.bam {temp_folder}/mapping_unpaired.bam 2>> {logfile}
    """, shell=True)
if method == "minimap2":
    call(f"""
    minimap2 -x sr --secondary=no -t 24  {temp_folder}/binset.idx {temp_folder}/fwd.fastq {temp_folder}/rev.fastq -a -o {temp_folder}/mapping.sam --MD  >> {logfile} 2>&1
    samtools view  --reference {temp_folder}/binset.fna -F0x900 -b -S -@{threads}  {temp_folder}/mapping.sam | samtools sort -@ 24 -o {temp_folder}/mapping_pairs.bam - >> {logfile} 2>&1
    minimap2 -x sr --secondary=no -t 24  {temp_folder}/binset.idx {temp_folder}/unp.fastq -a -o {temp_folder}/mapping.sam --MD  >> {logfile} 2>&1
    samtools view  --reference {temp_folder}/binset.fna -F0x900 -b -S -@{threads}  {temp_folder}/mapping.sam | samtools sort -@ 24 -o {temp_folder}/mapping_unpaired.bam - >> {logfile} 2>&1
    rm {temp_folder}/mapping.sam 2>> {logfile}
    samtools merge -f -t {threads} {temp_folder}/mapping.bam  {temp_folder}/mapping_pairs.bam  {temp_folder}/mapping_unpaired.bam >> {logfile}  2>&1
    #rm {temp_folder}/mapping_pairs.bam {temp_folder}/mapping_unpaired.bam 2>> {logfile}
    """, shell = True)


call(f"""
coverm filter -b {temp_folder}/mapping_pairs.bam -o {temp_folder}/mapping_filtered.bam --min-read-percent-identity {ani_cutoff} --min-read-aligned-length {min_len} --threads {threads} >> {logfile}  2>&1
#rm {temp_folder}/mapping.bam >> {logfile}  2>&1
""", shell=True)
title2log(f"extracting mapped reads from  {sample}", logfile)
call(f"samtools fastq -@ {threads} {temp_folder}/mapping_filtered.bam -1 {temp_folder}/matched_fwd.fastq -2 {temp_folder}/matched_rev.fastq -s {temp_folder}/matched_unp.fastq  2>> {logfile}", shell=True)

title2log("Done with the mappings", logfile)

title2log("Making tables", logfile)

shutil.copy(f"{temp_folder}/matched_fwd.fastq", f"{out_folder}/results/mappings/{map_name}/{sample}_fwd.fastq")
shutil.copy(f"{temp_folder}/matched_rev.fastq", f"{out_folder}/results/mappings/{map_name}/{sample}_rev.fastq")
shutil.copy(f"{temp_folder}/matched_unp.fastq", f"{out_folder}/results/mappings/{map_name}/{sample}_unp.fastq")

title2log("Cleaning up and moving the bins", logfile)
shutil.rmtree(temp_folder)
