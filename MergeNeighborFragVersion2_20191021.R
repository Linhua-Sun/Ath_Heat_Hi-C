## Aim: merge nearby fragments into segments
## Linhua Sun
## Mon Oct 21 16:32:26 CST 2019

library(ggpubr)
library(data.table)
library(rtracklayer)
library(magrittr)

setwd("/data1/linhua/QIANLAB/PROJECT/hic/fragments")

## Merge neighbouring fragments ----

RF<-fread("/data1/linhua/QIANLAB/PROJECT/hic/fragments/TAIR10_chr.tsv")

colnames(RF)[1:4] <- c("seqnames", "start", "end", "length")

MergeNeighborFrag <- function(data) {
    gg <- data
    GG <- gg ## Control gg is unchanged.
    i = 1
    while (i < nrow(gg)) {
        print(paste("Start ", i, sep = ""))
        if (GG[i]$length < 500) {
            j = 0
            while (GG[i]$length < 500) {
                j <- j + 1
                GG[i]$length <- GG[i]$length + gg[i + j]$length
            }
            GG[i]$end <- GG[i + j]$end
            for (k in 1:j) {
                GG[i + k]$length <- GG[i]$length
                GG[i + k]$start <- GG[i]$start
                GG[i + k]$end <- GG[i + j]$end
            }
            i <- i + k
        } else {
            i <- i + 1
        }
        print(paste("End ", i, sep = ""))
    }
    return(unique(GG[,1:4]))
}

Chr1_RF<-MergeNeighborFrag(data = RF[seqnames=="Chr1"])
Chr2_RF<-MergeNeighborFrag(data = RF[seqnames=="Chr2"])
Chr3_RF<-MergeNeighborFrag(data = RF[seqnames=="Chr3"])
Chr4_RF<-MergeNeighborFrag(data = RF[seqnames=="Chr4"])
Chr5_RF<-MergeNeighborFrag(data = RF[seqnames=="Chr5"])

rbind(Chr1_RF,
      Chr2_RF,
      Chr3_RF,
      Chr4_RF,
      Chr5_RF) -> New_RF

fwrite(New_RF,file = "AthMboI_MergedFragments.tsv",sep = "\t",col.names = F)

## Small bug in the tail of table ----
#       seqnames    start      end length   Group
#1:     Chr5 26971431 26972173    743 segment
#2:     Chr5 26972174 26973015    842 segment
#3:     Chr5 26973016 26973587    572 segment
#4:     Chr5 26973588 26974112    525 segment
#5:     Chr5 26974113 26975498   1386 segment
#6:     Chr5 26975499 26975502      4 segment


## Compare in silo fragment length distribution ----

New_RF<-fread("AthMboI_MergedFragments.tsv")
colnames(New_RF) <- c("seqnames","start","end","length","Group")

RF$Group<-"restriction fragment"
New_RF$Group<-"segment"

rbind(RF[, c("seqnames","start","end","length","Group")], New_RF[, c("seqnames","start","end","length","Group")])[,c("length","Group")][] -> temp

pdf(file = "Compare_RF_Merge.pdf",height = 6,width = 8)
gghistogram(
    temp[length <= 2000],
    x = "length",
    fill = "Group",
    palette = c("#00AFBB", "#E7B800")
) + border() + theme(aspect.ratio = 2 / (1 + sqrt(5))) 
dev.off()

## Import pericentromere regions (Chang Liu version) ----
peri<-fread("/data1/linhua/QIANLAB/PROJECT/hic/TAIR10_Pericentromere.txt")
peri$V2<-peri$V2*1000000
peri$V3<-peri$V3*1000000
colnames(peri)<-c("seqnames","start","end")
peri<-GRanges(peri)

## Remove pericentromere regions' intervales from the merged datasets ----

GRF<-GRanges(New_RF[,1:4])
GRF$ID<-seq(1,length(GRF))

overlapped<-findOverlapPairs(peri,GRF)

GRF[!(GRF$ID %in% overlapped@second$ID)] -> Clean_RF

as.data.table(Clean_RF)[,1:4] %>% fwrite("NoPR_AthMboI_MergedFragments.tsv",sep = "\t",row.names = F,col.names = F)

## Note: Mannul remove 'Chr5 26975499 26975502     4' merge into last term.

## Gererate IF Table ----

bed<-import("/data1/linhua/QIANLAB/PROJECT/hic/New/TAIR10_chr.bed")
strand(bed)<-"*"

Clean_RF<-fread("NoPR_AthMboI_MergedFragments.tsv")
colnames(Clean_RF)<-c("seqnames","start","end","length")
Clean_RF$name<-paste("RF_",seq(1,nrow(Clean_RF)),sep="")
Clean_RF<-GRanges(Clean_RF)

FF<-findOverlapPairs(bed,Clean_RF)

cbind(as.data.table(FF@first)[, c("seqnames", "start", "end", "name","width")],
      as.data.table(FF@second)[, c("seqnames", "start", "end", "name","width")]) %T>%
    fwrite(file = "RF_SF_Table.tsv",sep = "\t",col.names = F) %$%
    as.data.table(.)[, c(4, 9)][] %>% unique() %T>%
    fwrite(file = "RF_SF_ID_Table.tsv",sep = "\t",col.names = F) %>% as.data.table() -> KK

## Export into bed files----

II<-fread("RF_SF_Table.tsv")

newbed<-unique(II[,6:9])
colnames(newbed)<-c("seqnames","start","end","name")
export.bed(GRanges(newbed),con = "NoPR_SegmentsLevelFrags.bed")

OldRFBed<-unique(II[,1:4])
colnames(OldRFBed)<-c("seqnames","start","end","name")
export.bed(GRanges(OldRFBed),con = "RestrictionLevelFrags.bed")

## Color version ----

## Define a function to convert GRanges into bed12 with colors and name (Yes|No) 
ConvertGr2Bed <- function(BED, Prefix = "DefaultPrefix", RmName = T) {
    library(rtracklayer)
    NG <- BED
    NGBed <- sort(asBED(split(NG, NG$name)))
    
    if ((length(NGBed) %% 2)==0) {
        NGBed$itemRgb <- rep(c("lightgrey", "black"), length(NGBed) / 2)
    } else {
        NGBed$itemRgb <-
            c(rep(c("lightgrey", "black"), round(length(NGBed) / 2)), "lightgrey")
        
    }
    
    export.bed(sort(NGBed), con = paste(Prefix, "_Color.bed", sep = ""))
    if (RmName) {
        NGBed$name <- NULL
        export.bed(sort(NGBed),
                   con = paste(Prefix, "_NoName_Color.bed", sep = ""))
        
    }
}

ConvertGr2Bed(BED = GRanges(newbed),  Prefix ="Ath_MboI_NoPR_SegmentsLevelFrags" ,RmName = T)
ConvertGr2Bed(BED = GRanges(OldRFBed),Prefix ="Ath_MboI_RestrictionLevelFrags" ,RmName = T)

## ----

## RF[,1:4] # raw restriction fragments
## New_RF[,1:4] # merged nearby fragments into segments
## as.data.table(Clean_RF)[,1:4] # removed segments from pericentromere

