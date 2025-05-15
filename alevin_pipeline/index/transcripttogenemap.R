library(tidyverse)
library(data.table)
library(splitstackshape)
library(tidyr)
dmel_gtf <- fread(file="./dmel-all-UASStinger-r6.30.gtf", header=F, sep="\t", sep2 = ";")
dmel_gtf <- cSplit(dmel_gtf, "V9", ";")
dmel_gtf_fbid <- dmel_gtf[,c(9,11)]
dmel_gtf_fbid <- cSplit(dmel_gtf_fbid, splitCols = 1:2, " ")
dmel_gtf_fbid$V9_1_2 <- gsub("\"","",dmel_gtf_fbid$V9_1_2)
dmel_gtf_fbid$V9_3_2 <- gsub("\"","",dmel_gtf_fbid$V9_3_2)
dmel_gtf_fbid <- dmel_gtf_fbid[,c(2,4)]
#View(dmel_gtf_fbid)
dmel_gtf_fbid <- dmel_gtf_fbid %>% drop_na()
#View(dmel_gtf_fbid)
dmel_gtf_fbid <- unique(dmel_gtf_fbid[,c('V9_1_2','V9_3_2')])
#View(dmel_gtf_fbid)
dmel_gtf_fbid <- dmel_gtf_fbid[,c(2,1)]
#View(dmel_gtf_fbid)
write.table(dmel_gtf_fbid, file="dmel_all_UASStinger_r6.30_FBID_txp2gene.tsv",sep = "\t", quote=F, row.names = F, col.names = F )

dmel_gtf_sym <- dmel_gtf[,c(10,12)]
dmel_gtf_sym <- cSplit(dmel_gtf_sym, splitCols = 1:2, " ")
dmel_gtf_sym$V9_2_2 <- gsub("\"", "", dmel_gtf_sym$V9_2_2)
dmel_gtf_sym$V9_4_2 <- gsub("\"","", dmel_gtf_sym$V9_4_2)
dmel_gtf_sym <- dmel_gtf_sym[,c(2,4)]
dmel_gtf_sym <- dmel_gtf_sym %>% drop_na()
dmel_gtf_sym <- unique(dmel_gtf_sym[,c(1:2)])
dmel_gtf_sym <- dmel_gtf_sym[,c(2,1)]
write.table(dmel_gtf_sym, file="dmel_all_UASStinger_r6.30_Sym_txp2gene.tsv", sep="\t", quote=F, row.names = F, col.names = F)

dmel_gtf_fbtr_genesym <- dmel_gtf[,c(10,11)]
dmel_gtf_fbtr_genesym <- cSplit(dmel_gtf_fbtr_genesym, splitCols = 1:2, " ")
dmel_gtf_fbtr_genesym$V9_2_2 <- gsub("\"", "", dmel_gtf_fbtr_genesym$V9_2_2)
dmel_gtf_fbtr_genesym$V9_3_2 <- gsub("\"","", dmel_gtf_fbtr_genesym$V9_3_2)
dmel_gtf_fbtr_genesym <- dmel_gtf_fbtr_genesym[,c(2,4)]
dmel_gtf_fbtr_genesym <- dmel_gtf_fbtr_genesym %>% drop_na()
dmel_gtf_fbtr_genesym <- unique(dmel_gtf_fbtr_genesym[,c(1:2)])
dmel_gtf_fbtr_genesym <- dmel_gtf_fbtr_genesym[,c(2,1)]
write.table(dmel_gtf_fbtr_genesym, file="dmel_all_UASStinger_r6.30_FBTr_geneSym_txp2gene.tsv", sep="\t", quote=F, row.names = F, col.names = F)
