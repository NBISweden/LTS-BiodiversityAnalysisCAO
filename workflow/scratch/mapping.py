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
### This script is based in the mapping.py script from https://github.com/moritzbuck/0053_metasssnake2 commit 3ca68f087d7f0faa65c3400e14a9779cbb18b468


def title2log(title, logfile, llen = 90, also_stderr = True) :
    text_insert = "{title} started at : {time}".format(title = title, time = datetime.now())
    prefix = "="*floor((llen-2-len(text_insert))/2) + " "
    sufffix = " " + "="*ceil((llen-2-len(text_insert))/2)
    text = prefix + text_insert + sufffix
    with open(logfile, "a") as handle:
        handle.writelines("\n\n" + text + "\n\n")
    if also_stderr:
        print(text, file = stderr, flush = True)

def freetxt_line(text, logfile, llen = 90, also_stderr = True) :
    text_insert =  text
    prefix = "="*floor((llen-2-len(text_insert))/2) + " "
    sufffix = " " + "="*ceil((llen-2-len(text_insert))/2)
    text = prefix + text_insert + sufffix
    with open(logfile, "a") as handle:
        handle.writelines("\n\n" + text + "\n\n")

    if also_stderr:
        print(text, file = stderr, flush = True)


script , sample_file,reference, map_name, out_folder, temp_folder, mapper, ani_cutoff, min_len, threads, lib = sys.argv

sample_df = pandas.read_csv(sample_file, index_col=0)
samples = sample_df.to_dict(orient="index")

settings = {
    "reference" : reference,
    "out_folder" :  out_folder,
    "temp_folder" : temp_folder,
    "method" : mapper,
    "ani_cutoff" : ani_cutoff,
    "min_len" : min_len,
    "threads" : threads,
    "lib" :  lib,
    "fwds" : samples[wildcards.sample]["fwd_libs"],
    "revs" :  samples[wildcards.sample]["rev_libs"]
}

settings['temp_folder'] = mkdtemp(prefix = pjoin(settings['temp_folder'], "coi_mapping_"))
settings['logfile'] = pjoin(settings["out_folder"], 'logs', "run.log")

settings['fwds'] = settings['fwds'].split(";")
settings['revs'] = settings['revs'].split(";")

os.makedirs(os.path.dirname(settings['logfile']), exist_ok=True)

locals().update(settings)

call(f"conda env export > {out_folder}/logs/run_env.yaml", shell=True)
with open(pjoin(settings["out_folder"], 'logs', "run_settings.yaml"), "w") as handle:
    yaml.dump(settings)


os.makedirs(temp_folder, exist_ok=True)

title2log("reference to temp_folder", logfile)

shutil.copy(reference, pjoin(temp_folder, "ref.fna"))

title2log("indexing binset to temp_folder", logfile)
if method == "bwa-mem2":
    call(f"bwa-mem2 index {temp_folder}/ref.fna >> {logfile}  2>&1", shell=True)
if method == "bowtie2":
    call(f"bowtie2-build --threads {threads} {temp_folder}/ref.fna {temp_folder}/ref.fna >> {logfile}  2>&1", shell=True)
if method == "bbmap.sh":
    pass
if method == "minimap2":
       call(f"minimap2 -I 30G -t {threads} -d {temp_folder}/ref.idx {temp_folder}/ref.fna", shell=True)

freetxt_line(f"Starting to copy {len(fwds)} pairs of fastqs", logfile)

for ff, rr in zip(fwds,revs):
    title2log(f"copying {ff} and {rr} to temp_folder", logfile)
    call(f"""
    unpigz -kc {ff}> {temp_folder}/fwd.fastq 2>> {logfile}
    unpigz -kc {rr} > {temp_folder}/rev.fastq 2>> {logfile}
    """, shell=True)

call(f"fastp -y -l {min_len} -h /dev/null -j {out_folder}/fastp.json  --in1 {temp_folder}/fwd.fastq --in2 {temp_folder}/rev.fastq --out1 {temp_folder}/clean_fwd.fastq --out2 {temp_folder}/clean_rev.fastq -w {threads}  >> {logfile} 2>&1", shell= True)


freetxt_line("Starting mappings", logfile)
title2log("mapping {lib} to ref with {method}".format(lib = lib, method = method), logfile)
if method == "bwa-mem2":
    call(f"""
    bwa-mem2 mem -t {threads} {temp_folder}/ref.fna  -o {temp_folder}/mapping.sam {temp_folder}/clean_fwd.fastq {temp_folder}/clean_rev.fastq 2>> {logfile}
    samtools view  -b -S -@{threads}  {temp_folder}/mapping.sam | samtools sort -@ 24 -o {temp_folder}/mapping_pairs.bam - >> {logfile} 2>&1
    rm {temp_folder}/mapping.sam 2>> {logfile}
    """, shell=True)
if method == "bowtie2" :
    call(f"""
    bowtie2 -p {threads} -x  {temp_folder}/ref.fna --very-sensitive -S {temp_folder}/mapping.sam  -1 {temp_folder}/clean_fwd.fastq -2 {temp_folder}/clean_rev.fastq >> {logfile} 2>&1
    samtools view  -b -S -@{threads}  {temp_folder}/mapping.sam | samtools sort -@ 24 -o {temp_folder}/mapping_pairs.bam - >> {logfile} 2>&1
    bowtie2 -p {threads} -x  {temp_folder}/ref.fna --very-sensitive -S {temp_folder}/mapping.sam  -U {temp_folder}/unp.fastq >> {logfile} 2>&1
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
    rm {temp_folder}/mapping_pairs.bam {temp_folder}/mapping_unpaired.bam 2>> {logfile}
    """, shell = True)


call(f"""
coverm filter -b {temp_folder}/mapping_pairs.bam -o {temp_folder}/mapping_filtered.bam --min-read-percent-identity {ani_cutoff} --min-read-aligned-length {min_len} --threads {threads} >> {logfile}  2>&1
#rm {temp_folder}/mapping.bam >> {logfile}  2>&1
""", shell=True)
title2log("extracting mapped reads from  {lib}".format(lib = lib), logfile)
call(f"samtools fastq -@ {threads} {temp_folder}/mapping_filtered.bam -1 {temp_folder}/matched_fwd.fastq -2 {temp_folder}/matched_rev.fastq -s {temp_folder}/matched_unp.fastq  2>> {logfile}", shell=True)

title2log("Done with the mappings", logfile)

title2log("Making tables", logfile)

call(f"cp {temp_folder}/matched_fwd.fastq {out_folder}/{lib}_fwd.fastq", shell=True)
call(f"cp {temp_folder}/matched_rev.fastq {out_folder}/{lib}_rev.fastq", shell=True)
call(f"cp {temp_folder}/matched_unp.fastq {out_folder}/{lib}_unp.fastq", shell=True)

title2log("Cleaning up and moving the bins", logfile)
shutil.rmtree(temp_folder)

#shutil.rmtree(temp_folder)
