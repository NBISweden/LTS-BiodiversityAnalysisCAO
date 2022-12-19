#!/usr/bin/env Rscript
library("argparse")
library("dada2")
library("ShortRead")


parser = ArgumentParser(description= "This program takes a list of fastq files as input and return classification based on dada2's classifier")

parser$add_argument('--input', '-i', help= 'one fastq files', nargs = "+")
parser$add_argument('--output', '-o', help= 'tsv output', default = "/dev/stdout")
parser$add_argument('--cutoff', '-c', help= 'cutoff', type= 'double', default = 0.8)
parser$add_argument('--threads', '-t', help= 'number of threads', type= 'integer', default = 20)
parser$add_argument('--database', '-d', help= 'database')

xargs = parser$parse_args()

cutoff = xargs$cutoff
file = xargs$input 
ref = xargs$database
threads = xargs$threads
output = xargs$output
print(xargs)


print("Loading reads")
reads = readFastq(file)
seqs = as.character(reads@sread)
print("Classifying reads")

classification = assignTaxonomy(seqs = seqs, refFasta=ref, multithread=threads, minBoot = 0, outputBootstraps = TRUE) 

print("Parsing fixing")
ss = row.names(classification$tax)
col2 = sapply(ss, function(x) paste0(sub("__", ":", sub("()",")",paste(classification$tax[x,], classification$boot[x,]/100, ")", sep="("), fixed = TRUE)), collapse = ","))
col4 = sapply(ss, function(x) paste0(sub("__",":", classification$tax[x,classification$boot[x,]/100 > 0.05]), collapse=","))
col1 = as.character(reads@id)
col3 = rep(".", length(ss))

df = data.frame(id = col1, with_boot = col2, dir = col3, filtered = col4)
#classification = sapply(classification, function(x) paste0(x[!is.na(x)], collapse=","))

write.table(df, file = output, sep = "\t", quote = FALSE, col.names=FALSE)