library("ShortRead")
library("dada2")
library(data.table)
library(hues)
library("rjson")
library(ggplot2)

files = list.files("/home/moritz/temp/mappings/BOLDvsMIMEDNA/", full.names = TRUE)
files = files[grepl("MOSAIC.*.fastq", files)]

tt = fromJSON(file = "/home/moritz/temp/mappings/BOLDvsMIMEDNA/coverages.json")

classify = function( fastq ){
print(paste0("processing ", fastq))
reads = readFastq(fastq)
seqs = as.character(reads@sread)
classification = assignTaxonomy(seqs = seqs, refFasta="/home/moritz/dbs/proteins/coI/BOLD-clean_COI5P_Animalia-dada_reps.fasta", multithread=23)
row.names(classification) = as.character(reads@id)
classification = apply(classification,1, paste, sep=";", collapse=";")
return(table(classification))
}

classify_fasta = function( fasta ){
print(paste0("processing ", fasta))
reads = readFasta(fasta)
seqs = as.character(reads@sread)
classification = assignTaxonomy(seqs = seqs, refFasta="/home/moritz/dbs/proteins/coI/BOLD-clean_COI5P_Animalia-dada_reps.fasta", multithread=23)
row.names(classification) = as.character(reads@id)
classification = apply(classification,1, paste, sep=";", collapse=";")
return(classification)
}

assied_reads = classify_fasta("/home/moritz/temp/mappings/BOLDvsMIMEDNA/spades_assy/transcripts.fasta")

out2 = sapply(files, classify)
all_tax = names(table(unlist(sapply(out2, function(x) names(x)))))

tot_reads = tt$total_reads
tot_reads = t(as.data.frame(tot_reads, check.names = FALSE))

md = read.csv("~/projects/mosaic/M001_MIME/config_nd_metadata/mosaic_general_table_sample.csv")
md = md[!duplicated(md$sample_label),]
gps = read.csv("~/projects/mosaic/M001_MIME/config_nd_metadata/mosaic_gps.tsv", sep="\t")
row.names(md) = md$sample_label
row.names(md) = gsub("PS122_CN_DNA_","MOSAIC-MIME-DNA-",row.names(md))
md$nearest_gps = sapply(md$sampling_date, function(t) gps$date[which.min(sapply(gps$date, function(x) abs(as.numeric(as.Date(gps$date) - as.Date(t)))))])
row.names(gps) = gps$date
md = cbind(md, gps[md$nearest_gps,])

ab_table = sapply(out2, function(x) x[all_tax])
row.names(ab_table) = all_tax
colnames(ab_table) = gsub(".fastq", "", sapply(colnames(ab_table), basename))
ab_table[is.na(ab_table) ] = 0
ab_table = data.table(ab_table)
ab_table[, classif := all_tax]

ab_vals = melt(ab_table, id.vars = "classif", variable.name = "sample", value.name="nb_reads")
ab_vals[, class := sapply(strsplit(classif,";"),function(x) paste(x[1:5], collapse=";"))]
ab_vals = ab_vals[, .(nb_reads = sum(nb_reads)), by =  c('class', 'sample')]
ab_vals[, nb_reads := nb_reads/tot_reads[sample]]
ab_vals = cbind(ab_vals, md[as.character(ab_vals$sample),])
ab_vals2 = ab_vals[!is.na(ab_vals$sampling_date),]
bla = ab_vals[!grepl("NA", class),.(total_data=mean(nb_reads)) , by = class]
#bla = bla[!grepl("Lepidosa", class)]
#bla = bla[!grepl("Appendicularia", class)]

blabla = bla$total_data
names(blabla) = bla$class
top10 = names(sort(blabla, decreasing=TRUE)[1:10])
ab_vals2 = ab_vals2[class %in% top10 & nb_reads > 0]
ab_vals2$class = factor(ab_vals2$class, levels = top10)
ggplot(ab_vals2[class %in% top10], aes(x=sampling_date, y=nb_reads, fill=class))+geom_bar(stat = "identity")+scale_fill_manual(values=iwanthue(length(table(ab_vals$class))))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+facet_grid(~water.type, scale="free")+theme_minimal()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("sampling date")+ylab("relative abundance of reads")

#ggplot(ab_vals2, aes(x=lat, y=nb_reads, col=class,group=class))+geom_line()+scale_color_manual(values=iwanthue(length(table(ab_vals$class))))+facet_grid(water.type~class, scale="free")+ theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggplot(ab_vals, aes(x=lat, y=nb_reads, col=class,group=class))+geom_line()+scale_color_manual(values=iwanthue(length(table(ab_vals$class))))+ theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



"MOSAIC-MIME-DNA-1003",155829015.0
"MOSAIC-MIME-DNA-1203", 153183683.0
"MOSAIC-MIME-DNA-1303", 202057063.0
"MOSAIC-MIME-DNA-1503", 170755687.0
"MOSAIC-MIME-DNA-1603", 129311994.0
"MOSAIC-MIME-DNA-1703", 210584492.0
"MOSAIC-MIME-DNA-17", 134826675.0
"MOSAIC-MIME-DNA-25", 157717497.0
"MOSAIC-MIME-DNA-3003", 136471907.0
"MOSAIC-MIME-DNA-3103", 141927322.0
"MOSAIC-MIME-DNA-31", 165545551.0
"MOSAIC-MIME-DNA-3303", 175063943.0
"MOSAIC-MIME-DNA-603", 144570888.0
"MOSAIC-MIME-DNA-703", 173612306.0
"MOSAIC-MIME-DNA-903", 135598935.0
