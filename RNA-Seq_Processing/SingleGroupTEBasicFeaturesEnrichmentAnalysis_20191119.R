#!/usr/bin/env Rscript

## Aims ----
## Author: Linhua Sun
## Date: 2019-11-19 (Summary from previous studies)

# Single group of TE to compare
# Boxplot of TE length, GC contents.
# Barplot of TE classification (Class and Superfamily)
# Barplot of superfamily enrichment analysis + family level enrichment analysis
# Chromosome wide distribution of TEs (color based on groups)

## Set up input paramaters ----
library("optparse")
suppressPackageStartupMessages(library(crayon))
option_list = list(
    make_option(
        c("-f", "--file"),
        type = "character",
        default = NULL,
        help = "TE list file",
        metavar = "character"
    ),
    make_option(
        c("-o", "--out"),
        type = "character",
        default = "out.txt",
        help = "output file name [default= %default]",
        metavar = "character"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$file)){
    print_help(opt_parser)
    stop("At least one argument must be supplied: -f TEListFilePath.txt", call.=FALSE)
}

cat(yellow(paste("The input file is ", opt$file,  sep = "")),sep = "\n")
cat(yellow(paste("The output file prefix is ", opt$out, sep = "")),sep = "\n")


## Library R packages ----
library(pacman)
suppressWarnings(suppressPackageStartupMessages(
    p_load(
        ggpubr,
        ggplot2,
        ggsci,
        stringr,
        data.table,
        hashmap,
        crayon,
        GenomicRanges,
        gtrellis,
        circlize,
        magrittr,
        scales
    )
))
## Starting ----

#Example: opt$file<-"Heat_116_ddccm.txt"

# setwd("/data1/linhua/QIANLAB/PROJECT/te/BW/Typical_BW")
#opt<-NULL
#opt$file<-"Heat_Activated_TEs.txt"
#opt$out<-"Heat_Activated_TE_Finnal"

TE_Group<-str_remove_all(opt$file,".txt")

F1<-fread(opt$file,header = F,col.names = "ID")

TE_Annotation <- fread("/data1/linhua/QIANLAB/PROJECT/te/TE_Annotation-2019-11-19-PlusH1AndMNase.csv") # To be replaced by new annotation

left_join(F1,TE_Annotation,by="ID") %>% as.data.table() -> F1_ANN

F1_ANN[,"group":=TE_Group]

Total<-TE_Annotation[,"group":="Total"]

rbind(F1_ANN,Total) -> MM

#MM %>% group_by(group) %>% summarise(count=n())
MM[,c("Position", "seqnames","start","end","Ol_TEGs"):=NULL]

MM$group<-factor(MM$group,c("Total",TE_Group))

## calculate enrichment of blabla ----
EnrichmentRatiolength<-round(log2(mean(na.omit(MM[group!="Total"]$length))/mean(na.omit(MM[group=="Total"]$length))),digits = 3)
EnrichmentRatioGC<-round(log2(mean(na.omit(MM[group!="Total"]$score.GC))/mean(na.omit(MM[group=="Total"]$score.GC))),digits = 3)
EnrichmentRatioH3K9me2<-round(log2(mean(na.omit(MM[group!="Total"]$score.H3K9me2))/mean(na.omit(MM[group=="Total"]$score.H3K9me2))),digits = 3)

# round((mean(na.omit(MM[group!="Total"]$score.H1.1))/mean(na.omit(MM[group=="Total"]$score.H1.1))),digits = 3)
# round((mean(na.omit(MM[group!="Total"]$score.H1.2))/mean(na.omit(MM[group=="Total"]$score.H1.2))),digits = 3)
# round((mean(na.omit(MM[group!="Total"]$score.H1.3))/mean(na.omit(MM[group=="Total"]$score.H1.3))),digits = 3)
# round((mean(na.omit(MM[group!="Total"]$score.MNase))/mean(na.omit(MM[group=="Total"]$score.MNase))),digits = 3)

print(paste("Enrichment Ratio length is ",EnrichmentRatiolength,sep = ""))
print(paste("Enrichment Ratio GC content is ",EnrichmentRatioGC,sep = ""))
print(paste("Enrichment Ratio H3K9me2 is ",EnrichmentRatioH3K9me2,sep = ""))

fwrite(MM[group!="Total"],file = paste(TE_Group,"_AnnWithMultiFeatures.csv",sep=""))

rbind(
compare_means(data = MM[!is.na(score.H1.1)],formula =score.H1.1 ~ group ),
compare_means(data = MM[!is.na(score.H1.2)],formula =score.H1.2 ~ group ),
compare_means(data = MM[!is.na(score.H1.3)],formula =score.H1.3 ~ group ),
compare_means(data = MM[!is.na(score.MNase)],formula =score.MNase ~ group ),
compare_means(data = MM[!is.na(score.H3K9me2)],formula =score.H3K9me2 ~ group ),
compare_means(data = MM[!is.na(score.GC)],formula =score.GC ~ group ),
compare_means(data = MM[!is.na(length)],formula =length ~ group )) %>% fwrite(file = paste(TE_Group,"_CompareMeans.csv",sep=""))

fwrite(MM,file =paste(TE_Group,"_AnnWithMultiFeaturesPlusTotal.csv",sep=""))

## Set UP ploting parameters ----
source("/data1/linhua/QIANLAB/R_Functions/theme_linhua.R")

theme_set(theme_linhua(base_size = 18,base_family = "",legend = "none",border = F))
theme_update(panel.grid.minor.x = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.grid.major.x = element_blank(),
             panel.grid.major.y = element_blank(),
             axis.line.x = element_blank(),
             axis.title.x = element_blank(),
             #axis.text.x = element_blank(),
             #axis.ticks.x = element_blank(),
             panel.spacing.x = unit(0.5, "lines"),
             panel.spacing.y = unit(0.5, "lines"),
             strip.background = element_blank(),
             panel.border = element_blank(),
             axis.line.y = element_line(colour = "black"),
             aspect.ratio = 2
)
## GC contents ggplot ----

pdf(file = str_c(str_remove_all(opt$file,".txt"),"_GC.pdf",sep = ""),height = 10,width = 4)

ggplot(MM[!is.na(score.GC)], aes(x = group, y = score.GC , fill = group)) +
    stat_boxplot(geom = 'errorbar', width = 0.4) +
    geom_boxplot(outlier.colour = "black",width = 0.5,outlier.size = 1,notch = T,notchwidth = 0,outlier.shape = 1) +
    #stat_summary(fun.y = median,geom = "line",aes(group = 1),size = 1,color = "deepskyblue2") +
    #stat_summary(fun.y = median,geom = "point",shape = 21,size = 3,fill = "yellow",alpha = 0.8) +
    scale_fill_aaas(alpha = 1) +
    labs(x = "", y = "GC contents") +
    scale_y_continuous() +
    border() +
    theme(aspect.ratio = 2,
          axis.text.x = element_text(angle = 30, h = 1)) -> GC_Plot

YRange<-layer_scales(GC_Plot)$y$range$range
LLL<-as.list(as.data.frame(combn(levels(MM$group),m = 2),stringsAsFactors = F))
names(LLL)<-NULL

GC_Plot + stat_compare_means(
    comparisons = LLL,
    label.y = seq(mean(YRange), YRange[2], length.out = length(LLL)),
    hide.ns = F,
    color = "blue"
)
dev.off()

## (ChIP-on-chip) H1.1 H1.2 H1.3 and (NucMap) MNase-seq ----
pdf(file = str_c(str_remove_all(opt$file,".txt"),"_H1.1.pdf",sep = ""),height = 10,width = 4)

ggplot(MM[!is.na(score.H1.1)], aes(x = group, y = score.H1.1 , fill = group)) +
    stat_boxplot(geom = 'errorbar', width = 0.4) +
    geom_boxplot(outlier.colour = "black",width = 0.5,outlier.size = 1,notch = F,notchwidth = 0,outlier.shape = 1) +
    stat_summary(fun.y = median,geom = "line",aes(group = 1),size = 1,color = "deepskyblue2") +
    stat_summary(fun.y = median,geom = "point",shape = 21,size = 3,fill = "yellow",alpha = 0.8) +
    scale_fill_aaas(alpha = 1) +
    labs(x = "", y = "H1.1") +
    #scale_y_continuous(trans = log10_trans(),
    #                   breaks = trans_breaks("log10", function(x) 10^x),
    #                   labels = trans_format("log10", math_format(10^.x)))+
    border() +
    theme(aspect.ratio = 2,
          axis.text.x = element_text(angle = 30, h = 1)) -> H1.1_Plot

YRange<-layer_scales(H1.1_Plot)$y$range$range
LLL<-as.list(as.data.frame(combn(levels(MM$group),m = 2),stringsAsFactors = F))
names(LLL)<-NULL

H1.1_Plot + stat_compare_means(comparisons = LLL,label.y = seq(mean(YRange), YRange[2], length.out = length(LLL)),hide.ns = F,color = "blue")
dev.off()



pdf(file = str_c(str_remove_all(opt$file,".txt"),"_H1.2.pdf",sep = ""),height = 10,width = 4)

ggplot(MM[!is.na(score.H1.2)], aes(x = group, y = score.H1.2 , fill = group)) +
    stat_boxplot(geom = 'errorbar', width = 0.4) +
    geom_boxplot(outlier.colour = "black",width = 0.5,outlier.size = 1,notch = F,notchwidth = 0,outlier.shape = 1) +
    stat_summary(fun.y = median,geom = "line",aes(group = 1),size = 1,color = "deepskyblue2") +
    stat_summary(fun.y = median,geom = "point",shape = 21,size = 3,fill = "yellow",alpha = 0.8) +
    scale_fill_aaas(alpha = 1) +
    labs(x = "", y = "H1.2") +
    #scale_y_continuous(trans = log10_trans(),
    #                   breaks = trans_breaks("log10", function(x) 10^x),
    #                   labels = trans_format("log10", math_format(10^.x)))+
    border() +
    theme(aspect.ratio = 2,
          axis.text.x = element_text(angle = 30, h = 1)) -> H1.2_Plot

YRange<-layer_scales(H1.2_Plot)$y$range$range
LLL<-as.list(as.data.frame(combn(levels(MM$group),m = 2),stringsAsFactors = F))
names(LLL)<-NULL

H1.2_Plot + stat_compare_means(comparisons = LLL,label.y = seq(mean(YRange), YRange[2], length.out = length(LLL)),hide.ns = F,color = "blue")
dev.off()

pdf(file = str_c(str_remove_all(opt$file,".txt"),"_H1.3.pdf",sep = ""),height = 10,width = 4)

ggplot(MM[!is.na(score.H1.3)], aes(x = group, y = score.H1.3 , fill = group)) +
    stat_boxplot(geom = 'errorbar', width = 0.4) +
    geom_boxplot(outlier.colour = "black",width = 0.5,outlier.size = 1,notch = F,notchwidth = 0,outlier.shape = 1) +
    stat_summary(fun.y = median,geom = "line",aes(group = 1),size = 1,color = "deepskyblue2") +
    stat_summary(fun.y = median,geom = "point",shape = 21,size = 3,fill = "yellow",alpha = 0.8) +
    scale_fill_aaas(alpha = 1) +
    labs(x = "", y = "H1.3") +
    #scale_y_continuous(trans = log10_trans(),
    #                   breaks = trans_breaks("log10", function(x) 10^x),
    #                   labels = trans_format("log10", math_format(10^.x)))+
    border() +
    theme(aspect.ratio = 2,
          axis.text.x = element_text(angle = 30, h = 1)) -> H1.3_Plot

YRange<-layer_scales(H1.3_Plot)$y$range$range
LLL<-as.list(as.data.frame(combn(levels(MM$group),m = 2),stringsAsFactors = F))
names(LLL)<-NULL

H1.3_Plot + stat_compare_means(comparisons = LLL,label.y = seq(mean(YRange), YRange[2], length.out = length(LLL)),hide.ns = F,color = "blue")
dev.off()

pdf(file = str_c(str_remove_all(opt$file,".txt"),"_MNase.pdf",sep = ""),height = 10,width = 4)

ggplot(MM[!is.na(score.MNase)], aes(x = group, y = score.MNase , fill = group)) +
    stat_boxplot(geom = 'errorbar', width = 0.4) +
    geom_boxplot(outlier.colour = "black",width = 0.5,outlier.size = 1,notch = F,notchwidth = 0,outlier.shape = 1) +
    stat_summary(fun.y = median,geom = "line",aes(group = 1),size = 1,color = "deepskyblue2") +
    stat_summary(fun.y = median,geom = "point",shape = 21,size = 3,fill = "yellow",alpha = 0.8) +
    scale_fill_aaas(alpha = 1) +
    labs(x = "", y = "MNase") +
    #scale_y_continuous(trans = log10_trans(),
    #                   breaks = trans_breaks("log10", function(x) 10^x),
    #                   labels = trans_format("log10", math_format(10^.x)))+
    border() +
    theme(aspect.ratio = 2,
          axis.text.x = element_text(angle = 30, h = 1)) -> MNase_Plot

YRange<-layer_scales(MNase_Plot)$y$range$range
LLL<-as.list(as.data.frame(combn(levels(MM$group),m = 2),stringsAsFactors = F))
names(LLL)<-NULL

MNase_Plot + stat_compare_means(comparisons = LLL,label.y = seq(mean(YRange), YRange[2], length.out = length(LLL)),hide.ns = F,color = "blue")
dev.off()

## H3K9me2 ----

pdf(file = str_c(str_remove_all(opt$file,".txt"),"_H3K9me2.pdf",sep = ""),height = 10,width = 4)

ggplot(MM[!is.na(score.H3K9me2)], aes(x = group, y = score.H3K9me2 , fill = group)) +
    stat_boxplot(geom = 'errorbar', width = 0.4) +
    geom_boxplot(outlier.colour = "black",width = 0.5,outlier.size = 1,notch = F,notchwidth = 0,outlier.shape = 1) +
    stat_summary(fun.y = median,geom = "line",aes(group = 1),size = 1,color = "deepskyblue2") +
    stat_summary(fun.y = median,geom = "point",shape = 21,size = 3,fill = "yellow",alpha = 0.8) +
    scale_fill_aaas(alpha = 1) +
    labs(x = "", y = "H3K9me2") +
    scale_y_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)))+
    border() +
    theme(aspect.ratio = 2,
          axis.text.x = element_text(angle = 30, h = 1)) -> H3K9me2_Plot

YRange<-layer_scales(H3K9me2_Plot)$y$range$range
LLL<-as.list(as.data.frame(combn(levels(MM$group),m = 2),stringsAsFactors = F))
names(LLL)<-NULL

H3K9me2_Plot + stat_compare_means(comparisons = LLL,label.y = seq(mean(YRange), YRange[2], length.out = length(LLL)),hide.ns = F,color = "blue")
dev.off()
## Lengths ggplot ----

pdf(file = str_c(str_remove_all(opt$file,".txt"),"_Length_Distribution.pdf",sep = ""),height = 10,width = 4)
#MM$length_log10<-log10(MM$length)
ggplot(MM[!is.na(length)], aes(x = group, y = length , fill = group)) +
    stat_boxplot(geom = 'errorbar', width = 0.4) +
    geom_boxplot(outlier.colour = "black",width = 0.5,outlier.size = 1,notch = F,notchwidth = 0,outlier.shape = 1) +
    stat_summary(fun.y = median,geom = "line",aes(group = 1),size = 1,color = "deepskyblue2") +
    stat_summary(fun.y = median,geom = "point",shape = 21,size = 3,fill = "yellow",alpha = 0.8) +
    scale_fill_aaas(alpha = 1) +
    labs(x = "", y = "length (bp)") +
    scale_y_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)))+
    border() +
    theme(aspect.ratio = 2,axis.text.x = element_text(angle = 30, h = 1)) -> Length_Plot
YRange<-layer_scales(Length_Plot)$y$range$range

Length_Plot + stat_compare_means(comparisons = LLL,label.y = seq(mean(YRange), YRange[2], length.out = length(LLL)),hide.ns = F,color = "blue")
dev.off()

## TE Class Enrichment analysis ----
pdf(file = str_c(str_remove_all(opt$file,".txt"),"_Class_Enrichment.pdf",sep = ""),width = 8,height = 8)
MM %>% group_by(group,Class) %>% summarise(count=n()) %>% mutate(freq=count/sum(count)) %>%
    ggbarplot(
        x = "group",
        y = "freq",
        color = "black",
        fill = "Class",
        palette = "aaas",
        label = "count",
        width = 0.7,
        position = position_fill(),
        size = 0.5,
        lab.size = 4,
        lab.col = "white",
        lab.pos = "in"
    ) %>% ggpar(
        legend = "right",legend.title = "",
        xlab = "",x.text.angle = 45,
        ylab = "TE Class",
        font.x = c(18),
        font.y = c(18),
        font.tickslab = c(15),
        font.legend = c(15)
    ) +
    #border() +
    #scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))+ # https://stackoverflow.com/questions/13701347/force-the-origin-to-start-at-0
    #theme(aspect.ratio = 1/0.618) +
    theme(aspect.ratio = 2,axis.line.x =element_blank(),axis.ticks.x = element_blank())
dev.off()

## ----

Ath_Super<-fread("/data1/linhua/QIANLAB/PROJECT/Code/data/Ath_Superfamily_Sim.tsv")
Ath_Super$count<-NULL

pdf(file = str_c(str_remove_all(opt$file,".txt"),"_Fraction_Detailed_Superfamily.pdf",sep = ""),width = 8,height = 8)
left_join(MM,Ath_Super,by="Transposon_Super_Family") %>%
    group_by(group,TE_Super_Family) %>%
    summarise(count=n()) %>%
    group_by(group) %>%
    mutate(freq=count/sum(count)) %>%
    ungroup() %>%
    ggbarplot(
        x = "group",
        y = "freq",
        color = "black",
        fill = "TE_Super_Family",
        palette = "npg",alpha=0.8,
        width = 0.7,
        position = position_fill(),
        size = 0.5
    ) %>% ggpar(
        legend = "right",
        xlab = "",
        ylab = "TE",
        font.x = c(20),x.text.angle = 45,
        font.y = c(20),
        font.tickslab = c(18),
        font.legend = c(18)
    ) +
    #border() +
    #scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))+ # https://stackoverflow.com/questions/13701347/force-the-origin-to-start-at-0
    #theme(aspect.ratio = (1 + sqrt(5)) / 2) # golden ratio portrait
    theme(aspect.ratio =3,axis.line.x =element_blank(),axis.ticks.x = element_blank())
dev.off()

## Prepare for chromosome-wide plot ----

Ath_Chromatin<-fread("/data1/linhua/QIANLAB/PROJECT/Code/data/Ath_cytoband_Centromere-200K_gieStain.csv")
Ath_Chromatin$group<-0
Ath_Chromatin[gieStain=="gpos100"]$group<-1 ## centromere
Ath_Chromatin[gieStain=="gpos95"]$group<-2 ## PR
Ath_Chromatin[!(gieStain %in% c("gpos100","gpos95"))]$group<-3
col_fun = circlize::colorRamp2(seq(1,3),c("red","darkgrey","lightgrey"))
Chr_df = data.frame(c("Chr1","Chr2","Chr3", "Chr4", "Chr5"), c(1,1,1,1,1),c(30427671,19698289, 23459830, 18585056, 26975502))
#Chr_df34 <- data.frame(c("Chr3","Chr4"), c(1,1),c(23459830,18585056))
## Activated TE distribution in genome----
rbind(F1_ANN,Total) -> TT
TT %>% group_by(group) %>% summarise(count=n())
#TT[,c("Position", "seqnames","start","end","Ol_TEGs","score.H3K9me2"):=NULL]
TT$group<-factor(TT$group,c("Total",TE_Group))

TT[,c("seqnames","start","end","group")] -> GG
## TE distribution in genome ----
G1_Density <-
    genomicDensity(GG[GG$group == TE_Group],
                   window.size = 1000000,
                   overlap = F)
G1_Density$group<-TE_Group[1]


HH <- GRanges(TT)

pdf(str_c(str_remove_all(opt$file,".txt"),"_Genome_Distribution.pdf",sep = ""),height = 8,width = 14)

gtrellis_layout(
    data = Chr_df,
    n_track = 3,
    track_ylab = c("Length (log10)", "Density", "Chromatin"),
    asist_ticks = F,
    title = 'Distribution of activated transposons',
    track_ylim = c(c(min(log10(
        HH[HH$group != "Total"]$length
    )), max(log10(
        HH[HH$group != "Total"]$length
    ))), c(0, max(G1_Density$pct)), c(0, 1)),
    track_height = c(2, 2, 0.1),
    track_axis = c(T, T, F)
)

add_points_track(HH[HH$group == TE_Group],
                 log10(HH[HH$group == TE_Group]$length),
                 gp = gpar(col = pal_aaas(alpha = 1)(9)[2]),
                 size = unit(2, "mm"))

L1<-G1_Density[G1_Density$group==TE_Group,]

add_lines_track(
    L1,
    L1$pct,
    
    area = T,
    gp = gpar(col = "pink", fill = pal_aaas(alpha = 0.5)(9)[1])
)

add_points_track(L1, L1$pct,track = current_track(),size = unit(1, "mm"))
add_heatmap_track(Ath_Chromatin, Ath_Chromatin$group, fill = col_fun)

dev.off()
