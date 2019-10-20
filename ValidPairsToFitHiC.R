#!/usr/bin/env Rscript

## AIM:
## Shortcoming: Slow speed
## Linhua Sun: Sun Oct 20 17:07:34 2019

##############################################################################################################

library(optparse)
suppressWarnings(suppressMessages(library("crayon",character.only = TRUE,quietly = TRUE,warn.conflicts = FALSE)))
option_list = list(
    make_option(
        c("-b", "--bed"),
        type = "character",
        default = NULL,
        help = "segment level bed file name",
        metavar = "character"
    ),
    make_option(
        c("-i", "--interaction"),
        type = "character",
        default = NULL,
        help = "interaction file from validpairs name",
        metavar = "character"
    ),
    make_option(
        c("-o", "--out"),
        type = "character",
        default = "TEMP",
        help = "output dir name",
        metavar = "character"
    ),
    make_option(
        c("-c", "--chr"),
        type = "logical",
        default = F,
        help = "whether write data for each chromosome",
        metavar = "logical"
    ),
    make_option(
        c("-p", "--PR"),
        type = "logical",
        default = T,
        help = "whether remove countact counts from pericentromere",
        metavar = "logical"
    )
)

opt_parser = OptionParser(option_list = option_list)

opt = parse_args(opt_parser)

if (is.null(opt$bed)) {
    print_help(opt_parser)
    stop(red("At least 3 argument must be supplied: Rscriipt -b -i -o "), call.=FALSE)
}

cat(yellow(paste("The segment file is ", opt$bed,  sep = "")),sep = "\n")
cat(yellow(paste("The interaction file is ", opt$interaction, sep = "")),sep = "\n")
cat(yellow(paste("The output name is ", opt$out, sep = "")),sep = "\n")
cat(yellow(paste("Write data for each chromosome: ", as.character(opt$chr), sep = "")),sep = "\n")
cat(yellow(paste("Remove countact counts from pericentromere: ", as.character(opt$PR), sep = "")),sep = "\n")


## generate interactions by `cat ../../hic/viewpoint/sample_allValidPairs|cut -f9-10|uniq -c|awk '{OFS="\t"; print $1,$2,$3 }' > sample_allValidPairs.counts`
##############################################################################################################

# See whether these packages exist on comp. If not, install.
package_list <- c("data.table",
                  "dplyr",
                  "rtracklayer",
                  "magrittr",
                  "gdata",
                  "hashmap",
                  "crayon")
cat(yellow("Import R packages!"),sep = "\n")
for(p in package_list) {
    if (!suppressWarnings(suppressMessages(require(
        p,
        character.only = TRUE,
        quietly = TRUE,
        warn.conflicts = FALSE
    )))) {
        install.packages(p, repos = "http://cran.r-project.org")
        suppressWarnings(suppressMessages(library(
            p,
            character.only = TRUE,
            quietly = TRUE,
            warn.conflicts = FALSE
        )))
    }
}

## use pacman as alternative to import a series R packages
## Main ----

#opt<-NULL
#opt$bed<-"New_TAIR10_Merged.bed"
#opt$interaction<-"test.counts"
#opt$out<-"Heat"
#Chr<-F

BED<-opt$bed
INTERACTION<-opt$interaction
OutName<-opt$out
Chr<-opt$chr
PR<-opt$PR

#setwd("/data1/linhua/QIANLAB/PROJECT/hic/New")
#Args<-NULL
#Args[1]<-"New_TAIR10_Merged.bed"
#Args[2]<-"Heat_RF_Counts.txt" ## How the file generated?
#Args[3]<-"Heat"

## Import bed of merged nearby fragments ----
bed<-as.data.table(import(BED))
bed1<-bed
bed2<-bed
colnames(bed1)<-paste(colnames(bed1),"1",sep = "")
colnames(bed2)<-paste(colnames(bed2),"2",sep = "")
bed1$mid1<-round((bed1$start1+bed1$end1)/2)
bed2$mid2<-round((bed2$start2+bed2$end2)/2)

## Import interactions between fragments ----

cat(yellow("Import merged bed file!"),sep = "\n")

txt<-fread(INTERACTION,header = F)

cat(yellow("Import interaction data from validpairs!"),sep = "\n")

## Import hash table from ID_table.tsv [ fragment ~ segment ] ----
KK<-fread("/data1/linhua/QIANLAB/PROJECT/hic/New/ID_Table.tsv",header = F,col.names = c("old","new")) ## Need to change it based on env.
HH<-hashmap(KK$old,KK$new)

## Sub old raw fragments by new segments ----

txt$V4<-HH[[txt$V2]]
txt$V5<-HH[[txt$V3]]

txt[,c("V1","V4","V5")] %>% 
    group_by(V4,V5) %>%
    summarise(SS=sum(V1)) %$%
    as.data.table(.)[,c("SS","V4","V5")] -> TXT

cat(yellow("Sub with data from ID_Table.tsv using hash!"),sep = "\n")

OldNames<-colnames(TXT)
NewNames<-c("count","name1","name2")
setnames(TXT,OldNames,NewNames)

## Join segment level contact counts and bed1 and bed2 ----

OutFileName <- paste(OutName,"_Int",".tsv",sep="")

inner_join(TXT, bed1, by = "name1") %>% 
    inner_join(bed2, by = "name2") %>%
    as.data.table() %$%
    .[, c("seqnames1","mid1","seqnames2","mid2","count")][order(seqnames1,mid1,seqnames2,mid2)] %T>%
    fwrite(file = OutFileName,sep="\t",col.names = F) %>%
    as.data.table() -> UU

system( paste("gzip ",OutFileName ,sep = "") )

cat(yellow("Converting!"),sep = "\n")

if (Chr) {
    cat(yellow("Write countact counts at segment level for each chromosome!"),sep = "\n")
    UU[seqnames1=="Chr1" & seqnames2=="Chr1"] %>% fwrite(file = paste(OutName,"_INT_Chr1",".tsv",sep=""),sep = "\t",col.names = F,nThread = 10)
    UU[seqnames1=="Chr2" & seqnames2=="Chr2"] %>% fwrite(file = paste(OutName,"_INT_Chr2",".tsv",sep=""),sep = "\t",col.names = F,nThread = 10)
    UU[seqnames1=="Chr3" & seqnames2=="Chr3"] %>% fwrite(file = paste(OutName,"_INT_Chr3",".tsv",sep=""),sep = "\t",col.names = F,nThread = 10)
    UU[seqnames1=="Chr4" & seqnames2=="Chr4"] %>% fwrite(file = paste(OutName,"_INT_Chr4",".tsv",sep=""),sep = "\t",col.names = F,nThread = 10)
    UU[seqnames1=="Chr5" & seqnames2=="Chr5"] %>% fwrite(file = paste(OutName,"_INT_Chr5",".tsv",sep=""),sep = "\t",col.names = F,nThread = 10)
    cat(yellow("Fast write interaction data!"),sep = "\n")
}

## Fragments data----

UU %>% group_by(seqnames1, mid1) %>% summarise(s1 = sum(count)) %>% as.data.table() -> S1
UU %>% group_by(seqnames2, mid2) %>% summarise(s1 = sum(count)) %>% as.data.table() -> S2

colnames(S1)<-c("seqnames","mid","count1")
colnames(S2)<-c("seqnames","mid","count2")

setkeyv(S1,c("seqnames","mid"))
setkeyv(S2,c("seqnames","mid"))

merge(S1,S2,all=T)[] -> Sc
Sc<-NAToUnknown(Sc,unknown = 0)


OutFileName <- paste(OutName,"_Frag",".tsv",sep="")

Sc[,count:=count1+count2][,column2:=0][,column5:=0][,c("seqnames","column2","mid","count","column5")][] %T>%
    fwrite(file = OutFileName,sep = "\t",col.names = F) %>% 
    as.data.table() -> SS

system( paste("gzip ",OutFileName ,sep = "") )

if (Chr) {
    SS[seqnames=="Chr1"] %>% fwrite(file = paste(OutName,"_Frag_Chr1",".tsv",sep=""),sep = "\t",col.names = F)
    SS[seqnames=="Chr2"] %>% fwrite(file = paste(OutName,"_Frag_Chr2",".tsv",sep=""),sep = "\t",col.names = F)
    SS[seqnames=="Chr3"] %>% fwrite(file = paste(OutName,"_Frag_Chr3",".tsv",sep=""),sep = "\t",col.names = F)
    SS[seqnames=="Chr4"] %>% fwrite(file = paste(OutName,"_Frag_Chr4",".tsv",sep=""),sep = "\t",col.names = F)
    SS[seqnames=="Chr5"] %>% fwrite(file = paste(OutName,"_Frag_Chr5",".tsv",sep=""),sep = "\t",col.names = F)
    #SS[seqnames=="Chr1"][,count:=1][] %>% fwrite(file = paste(OutName,"_Frag_Chr1",".tsv",sep=""),sep = "\t",col.names = F)
    cat(yellow("Fast write frqgments data!"),sep = "\n")
}


## Filter chromatin interactions and fragments (segments) based on location: pericentromere ----
if (PR) {
    INT <- UU
    BIN <- SS
    rm(UU,SS)
    options(scipen = 999)
    
    INT[!(
        mid1 %in% c(14800502, 3975226, 12799638, 3475472, 12174946) |
            mid2 %in% c(14800502, 3975226, 12799638, 3475472, 12174946)
    )][order(seqnames1,mid1,seqnames2,mid2)] -> Clean_INT
    
    Clean_INT <- Clean_INT[, lapply(.SD, format, scientific = FALSE)]
    
    OutFileName <- paste(OutName,"_NoPRs_INT.tsv",sep = "")
    fwrite(
        Clean_INT,
        file = OutFileName,
        sep = "\t",
        col.names = F,
        nThread = 10
    )
    system( paste("gzip ",OutFileName ,sep = "") )
    
    #_________________________________________________________________________________________________________
    BIN[!(mid %in% c(14800502, 3975226, 12799638, 3475472, 12174946))][order(seqnames,mid)] -> NoBIN
    
    NoBIN <- NoBIN[, lapply(.SD, format, scientific = FALSE)]
    
    OutFileName <- paste(OutName,"_NoPRs_FRAG.tsv",sep = "")
    fwrite(
        NoBIN,
        file = OutFileName,
        sep = "\t",
        col.names = F,
        nThread = 10
    )
    system( paste("gzip ",OutFileName ,sep = "") )
}
