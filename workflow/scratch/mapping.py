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


#script , mapping_name, config_file, root_folder, out_folder, threads = sys.argv
test_lib_folder = "/home/moritz/data/raw_data/mosaic/UL-3163/220408_A00181_0471_AHF7CCDSX3/220408_A00181_0471_AHF7CCDSX3/Sample_UL-3163-SO21-DNA-007/"
settings = {
    "reference" : "/home/moritz/projects/mosaic/M002_EFICA/data/BOLD-clean_COI5P_Animalia-dada.fasta",
    "out_folder" :  "/home/moritz/temp/mapping/test/",
    "temp_folder" : "/home/moritz/temp/",
    "method" : "bwa-mem2",
    "ani_cutoff" : 0.95,
    "min_len" : 80,
    "threads" : 20,
    "fwds" : ";".join([pjoin(test_lib_folder, f) for f in os.listdir(test_lib_folder) if "_R1_" in f and f.endswith(".fastq.gz")]),
    "revs" : ";".join([pjoin(test_lib_folder, f) for f in os.listdir(test_lib_folder) if "_R2_" in f and f.endswith(".fastq.gz")])
}

settings['temp_folder'] = mkdtemp(prefix = pjoin(settings['temp_folder'], "coi_mapping_"))
settings['logfile'] = pjoin(settings["out_folder"], 'logs', "run.log")

settings['fwds'] = settings['fwds'].split(";")
settings['revs'] = settings['revs'].split(";")

os.makedirs(os.path.dirname(settings['logfile']), exist_ok=True)

call(f"conda env export > {out_folder}/logs/run_env.yaml", shell=True)
with open(pjoin(settings["out_folder"], 'logs', "run_settings.yaml"), "w") as handle:
    yaml.dump(settings)

locals().update(settings)

os.makedirs(temp_folder, exist_ok=True)

title2log("reference to temp_folder", logfile)

shutil.copy(reference, pjoin(temp_folder, "ref.fna"))

title2log("indexing binset to temp_folder", logfile)
if method == "bwa-mem2":
    call(f"bwa-mem2 index {temp_folder}/ref.fna >> {logfile}  2>&1", shell=True)
    pass
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
    call("""
    bowtie2 -p {threads} -x  {temp_folder}/ref.fna --very-sensitive -S {temp_folder}/mapping.sam  -1 {temp_folder}/fwd.fastq -2 {temp_folder}/rev.fastq >> {logfile} 2>&1
    samtools view  -b -S -@{threads}  {temp_folder}/mapping.sam | samtools sort -@ 24 -o {temp_folder}/mapping_pairs.bam - >> {logfile} 2>&1
    bowtie2 -p {threads} -x  {temp_folder}/ref.fna --very-sensitive -S {temp_folder}/mapping.sam  -U {temp_folder}/unp.fastq >> {logfile} 2>&1
    samtools view -b -S -@{threads}  {temp_folder}/mapping.sam | samtools sort -@ 24 -o {temp_folder}/mapping_unpaired.bam - >> {logfile} 2>&1
    rm {temp_folder}/mapping.sam 2>> {logfile}
    samtools merge -f -t {threads} {temp_folder}/mapping.bam  {temp_folder}/mapping_pairs.bam  {temp_folder}/mapping_unpaired.bam >> {logfile}  2>&1
    rm {temp_folder}/mapping_pairs.bam {temp_folder}/mapping_unpaired.bam 2>> {logfile}
    """.format(root_folder = root_folder, lname = lib, temp=temp_folder, threads = threads, out_folder = out_folder), shell=True)
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


    call("""
    coverm filter -b {temp_folder}/mapping.bam -o {temp_folder}/mapping_filtered.bam --min-read-percent-identity {ani} --min-read-aligned-length {min_len} --threads {threads} >> {logfile}  2>&1
    rm {temp_folder}/mapping.bam >> {logfile}  2>&1
    coverm contig  --bam-files {temp_folder}/mapping_filtered.bam  --methods count  --threads {threads} > {temp_folder}/mapping.tsv 2>> {logfile}
    """.format(temp=temp_folder, threads = threads, ani = ani,  out_folder = out_folder,  min_len = min_len), shell=True)
    if keep_mapped:
        title2log("extracting mapped reads from  {lib}".format(lib = lib), logfile)
        call(f"samtools fastq -@ {threads} {temp_folder}/mapping_filtered.bam -o {temp_folder}/{lib}.fastq  2>> {logfile}", shell=True)
    if keep_unmapped:
        title2log("extracting mapped reads from  {lib}".format(lib = lib), logfile)
        call(f"samtools fastq -f 4 -@ {threads} {temp_folder}/mapping_filtered.bam -o {temp_folder}/{lib}.fastq  2>> {logfile}", shell=True)
    call("""
    rm {temp_folder}/mapping_filtered.bam
    cat {temp_folder}/fwd.fastq {temp_folder}/rev.fastq {temp_folder}/unp.fastq | wc -l > {temp_folder}/total_reads_x4.txt
    """.format(temp=temp_folder, threads = threads, ani = ani,  out_folder = out_folder,  min_len = min_len), shell=True)

        call("""
    rm {temp_folder}/fwd.fastq
    rm {temp_folder}/rev.fastq
    rm {temp_folder}/unp.fastq
    """.format(temp = temp_folder), shell = True)

        with open(pjoin(temp_folder, "mapping.tsv")) as handle:
            handle.readline()
            coverages[lib] = {}
            for l in handle:
                ll = l.strip().split()
                coverages[lib][ll[0]] = float(ll[1])
        with open(pjoin(temp_folder, "total_reads_x4.txt")) as handle:
            total_reads[lib] = int(handle.readline().strip())/4

        with open(f"{temp_folder}/coverages.json", "w") as handle:
            json.dump({"coverages" : coverages, "total_reads" : total_reads}, handle)



title2log("Done with the mappings", logfile)

title2log("Making tables", logfile)

if not is_rna :
    for k in total_reads:
        coverages[k]['unmapped'] = total_reads[k] - sum(coverages[k].values())

with open(pjoin(temp_folder, "total_reads_to_map.csv"), "w") as handle:
    handle.writelines(["library,total_reads\n"] + [f"{k},{v}\n" for k,v in total_reads.items()])

coverages = pandas.DataFrame.from_dict(coverages)
coverages.to_csv(pjoin(temp_folder, "contigs_mapped_reads.csv"), index_label = "contig_name")
normed_cov = coverages/coverages.sum()
normed_cov.to_csv(pjoin(temp_folder, "contigs_relative_abundance.csv"), index_label = "contig_name")

if is_rna :
    gene_stats  = pandas.read_csv(pjoin(root_folder, "binsets", binset, alternate_root + "_basics.csv" ), index_col=0)
    lens = [len(gene_stats.loc[g, 'representative_nucls']) for g in normed_cov.index]
    rpk = 10_000*(coverages.transpose()/lens).transpose()
    norm_facts = rpk.sum()/1_000_000
    tpm = rpk/norm_facts
    tpm.to_csv(pjoin(temp_folder, "tpm.csv"), index_label = "derep_gene")

    for l in ['root_eggNOG', 'COG_category', 'symbol', 'KO']:
        relabs = tpm.copy()
        relabs[l] = [t for t in gene_stats.loc[tpm.index,l] ]
        if l in ['KO', 'COG_category']:
            if l == "KO":
                relabs[l] = [ "" if v != v else v.split(",") for v in relabs[l]]
            else :
                relabs[l] = [ "" if v != v else  list(v) for v in relabs[l]]
            multiples = [i for i,n in zip(relabs.index,relabs[l]) if n == n and len(n) >1 ]
            lines = []
            for gc in tqdm(multiples):
                line = relabs.loc[gc]
                values = line[l]
                del line[l]
                line = list(line/len(values))
                lines += [ line + [values[j]] for j in range(len(values))]
            relabs = relabs.loc[[i for i,n in zip(relabs.index,relabs[l]) if len(n) == 1 or n == ""]]
            relabs[l] = ["NA" if len(v) == 0 else v[0] for v in relabs[l]]
            lines =  pandas.DataFrame(lines)
            lines.columns = relabs.columns
            relabs = pandas.concat([relabs,lines])

        relabs = relabs.groupby(l).sum()
        relabs.to_csv(pjoin(temp_folder, l + "_relative_abundance.csv"), index_label = l)


else :
    coverages['bin'] = ["_".join(c.split("_")[:-1]) if c != "unmapped"  else "unmapped" for c in coverages.index]
    coverages_by_bin = coverages.groupby("bin").sum()

    relab_by_bin = (coverages_by_bin/coverages_by_bin.sum())
    relab_by_bin.to_csv(pjoin(temp_folder, "bins_relative_abundance.csv"), index_label = "bin_name")

    levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    bin_stats  = csv2dict(pjoin(root_folder, "binsets", binset, binset + "_basics.csv" ))


    taxas = [bin_stats.get(bin_, {taxfield : 'unmapped'})[taxfield] for bin_ in relab_by_bin.index]
    for l in range(7):
        relabs = relab_by_bin.copy()
        relabs[levels[l]] = [";".join(t.split(";")[0:(l+1)]) for t in taxas]
        relabs = relabs.groupby(levels[l]).sum()
        relabs.to_csv(pjoin(temp_folder, levels[l] + "_relative_abundance.csv"), index_label = levels[l])

    relabs = relab_by_bin.copy()
    relabs["mOTU"] = [bin_stats.get(bin_, {'mOTU' : 'unmapped'})['mOTU'] for bin_ in relabs.index]
    relabs = relabs.groupby("mOTU").sum()
    relabs.to_csv(pjoin(temp_folder, "mOTU_relative_abundance.csv"), index_label = levels[l])


call(f"cp {temp_folder}/*.csv {root_folder}/mappings/{mapping_name}/", shell=True)
call(f"mkdir -p {root_folder}/mappings/{mapping_name}/mapped_reads/; cp {temp_folder}/*.fastq {root_folder}/mappings/{mapping_name}/mapped_reads/", shell=True)

title2log("Cleaning up and moving the bins", logfile)


#shutil.rmtree(temp_folder)
