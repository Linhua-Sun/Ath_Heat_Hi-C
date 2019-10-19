#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(gdata))
suppressPackageStartupMessages(library(hashmap))
suppressPackageStartupMessages(library(crayon))
cat(yellow("Import R packages!"),sep = "\n")
Args<-commandArgs(trailingOnly = T)

if (length(Args)==0) {
    stop("At least 1 argument must be supplied: Rscriipt script.R bed txt", call.=FALSE)
}

## Args[1]<-"TAIR10_chr.bed"
## Args[2]<-"test.txt"
## Args[3]<-"Heat"
bed<-as.data.table(import(Args[1]))
bed1<-bed
bed2<-bed
colnames(bed1)<-paste(colnames(bed1),"1",sep = "")
colnames(bed2)<-paste(colnames(bed2),"2",sep = "")
bed1$mid1<-round((bed1$start1+bed1$end1)/2)
bed2$mid2<-round((bed2$start2+bed2$end2)/2)

cat(yellow("Import merged bed file!"),sep = "\n")

txt<-fread(Args[2])

cat(yellow("Import interaction data from validpairs!"),sep = "\n")

KK<-fread("ID_Table.tsv",header = F,col.names = c("old","new")) ## Need to change it based on env.
HH<-hashmap(KK$old,KK$new)
txt$V4<-HH[[txt$V2]]
txt$V5<-HH[[txt$V3]]

txt[,c("V1","V4","V5")] %>% 
    group_by(V4,V5) %>%
    summarise(SS=sum(V1)) %$%
    as.data.table(.)[,c("SS","V4","V5")] -> txt

cat(yellow("Sub with data from ID_Table.tsv using hash!"),sep = "\n")

colnames(txt)<-c("count","name1","name2")

inner_join(txt, bed1, by = "name1") %>% 
    inner_join(bed2, by = "name2") %>%
    as.data.table() %$%
    .[, c("seqnames1","mid1","seqnames2","mid2","count")][order(seqnames1,mid1,seqnames2,mid2)] %T>%
    fwrite(file = paste(Args[3],"_Int",".tsv",sep=""),sep="\t",col.names = F) %>%
    as.data.table() -> UU

cat(yellow("Converting!"),sep = "\n")

colnames(UU)<-c("V1","V2","V3","V4","V5")

#UU<-fread("zcat /data1/linhua/QIANLAB/PROJECT/hic/Test_New/txt.txt.gz")

UU[V1=="Chr1" & V3=="Chr1"] %>% fwrite(file = paste(Args[3],"_INT_Chr1",".tsv",sep=""),sep = "\t",col.names = F,nThread = 10)
UU[V1=="Chr2" & V3=="Chr2"] %>% fwrite(file = paste(Args[3],"_INT_Chr2",".tsv",sep=""),sep = "\t",col.names = F,nThread = 10)
UU[V1=="Chr3" & V3=="Chr3"] %>% fwrite(file = paste(Args[3],"_INT_Chr3",".tsv",sep=""),sep = "\t",col.names = F,nThread = 10)
UU[V1=="Chr4" & V3=="Chr4"] %>% fwrite(file = paste(Args[3],"_INT_Chr4",".tsv",sep=""),sep = "\t",col.names = F,nThread = 10)
UU[V1=="Chr5" & V3=="Chr5"] %>% fwrite(file = paste(Args[3],"_INT_Chr5",".tsv",sep=""),sep = "\t",col.names = F,nThread = 10)

cat(yellow("Fast write interaction data!"),sep = "\n")


UU %>% group_by(V1, V2) %>% summarise(s1 = sum(V5)) %>% as.data.table() -> S1
UU %>% group_by(V3, V4) %>% summarise(s1 = sum(V5)) %>% as.data.table() -> S2

colnames(S1)<-c("seqnames","mid","count1")
colnames(S2)<-c("seqnames","mid","count2")
setkeyv(S1,c("seqnames","mid"))
setkeyv(S2,c("seqnames","mid"))

merge(S1,S2,all=T)[] -> Sc
Sc<-NAToUnknown(Sc,unknown = 0)
Sc[,count:=count1+count2][,column2:=0][,column5:=0][,c("seqnames","column2","mid","count","column5")][] %T>%
    fwrite(file = paste(Args[3],"_Frag",".tsv",sep=""),sep = "\t",col.names = F) %>% 
    as.data.table() -> SS

SS[seqnames=="Chr1"] %>% fwrite(file = paste(Args[3],"_Frag_Chr1",".tsv",sep=""),sep = "\t",col.names = F)
SS[seqnames=="Chr2"] %>% fwrite(file = paste(Args[3],"_Frag_Chr2",".tsv",sep=""),sep = "\t",col.names = F)
SS[seqnames=="Chr3"] %>% fwrite(file = paste(Args[3],"_Frag_Chr3",".tsv",sep=""),sep = "\t",col.names = F)
SS[seqnames=="Chr4"] %>% fwrite(file = paste(Args[3],"_Frag_Chr4",".tsv",sep=""),sep = "\t",col.names = F)
SS[seqnames=="Chr5"] %>% fwrite(file = paste(Args[3],"_Frag_Chr5",".tsv",sep=""),sep = "\t",col.names = F)

#SS[seqnames=="Chr1"][,count:=1][] %>% fwrite(file = paste(Args[3],"_Frag_Chr1",".tsv",sep=""),sep = "\t",col.names = F)
cat(yellow("Fast write frqgments data!"),sep = "\n")
