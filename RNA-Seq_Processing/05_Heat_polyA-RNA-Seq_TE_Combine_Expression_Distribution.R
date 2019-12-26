## ----
library(gtrellis)
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(gtrellis))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(rtracklayer))
setwd("/data1/linhua/QIANLAB/PROJECT/hic/HOMER/HiC_RNA_Seq_Heat_CSV")
## Import data ----
DLR<-fread("HeatTagDir_DLR_BIN10K_limmted_0.05.csv")
colnames(DLR)[1]<-"seqnames"

ICF<-fread("HeatTagDir_ICF_BIN10K_limmted_0.2.csv")
colnames(ICF)[1]<-"seqnames"

Chr_df = data.frame(
    seqnames=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),
    start=c(1, 1, 1, 1, 1),
    end=c(30427671, 19698289, 23459830, 18585056, 26975502)
)

Chr_df = data.frame(
    seqnames=c("Chr1", "Chr2"),
    start=c(1, 1),
    end=c(30427671, 19698289)
)

Chr_df = data.frame(
    seqnames=c("Chr3", "Chr4"),
    start=c(1, 1),
    end=c(23459830, 18585056)
)

#Chr_df<-fread("Ath_Cytoband_KEE_Version.csv")

## ----
TT<-fread("/data1/linhua/QIANLAB/PROJECT/te/HeatTEReanalysis2019-06-03/Updated_Heat_Stress_polyA-RNA-Seq_DEG_TEs_Ann_Results.csv")
UP<-TT[Heat_vs_Control=="UP"]
UP$log10<-log10(UP$Heat.mean.cpm - UP$Control.mean.cpm)
GRUP<-GRanges(UP)
DGRUP<-as.data.table(GRUP)
DD<-DGRUP[,c("seqnames","start","end","log10","Chromatin")]
colnames(DD)<-c("chr","start","end","value","color")
DD[color=="chromosome"]$color<-"a"
DD[color!="a"]$color<-"b"

UPTE<-DD

A<-UPTE[color=="b"] ## heterochromatin
B<-UPTE[color=="a"] ## euchromatin

Ath_Chromatin<-fread("/data1/linhua/QIANLAB/PROJECT/hic/HOMER/BEDG/Ath_cytoband_Centromere-200K_gieStain.csv")
Ath_Chromatin$group<-0
Ath_Chromatin[gieStain=="acen"]$group<-1
Ath_Chromatin[gieStain=="gneg"]$group<-2
Ath_Chromatin[gieStain=="gpos100"]$group<-3
Ath_Chromatin[gieStain=="gpos95"]$group<-4
Ath_Chromatin[gieStain=="stalk"]$group<-5
col_fun = circlize::colorRamp2(seq(1,5),c("red","white","black","grey","darkblue"))

mypal<-pal_npg(alpha = 1)(9)

## ----

pdf(file = "UP_gtrellis_Chr34_ICF_DLR.pdf",height = 6,width = 10)
#range(UPTE[chr %in% c("Chr3","Chr4")]$value)
gtrellis_layout(
    data = Chr_df,
    gap = 0,
    track_ylim = c(c(0,1),c(-0.5,2.2), c(-0.05, 0.05), c(-0.2, 0.2)),
    asist_ticks = T,
    add_name_track = T,
    name_fontsize = 12,
    #legend = lgd,
    name_track_fill = "white",
    track_ylab = c("","TE", "DLR", "ICF"),
    axis_label_fontsize = 8,
    lab_fontsize = 12,
    remove_chr_prefix = F,
    n_track = 4,
    track_height = c(0.2,1,1, 1),
    track_axis = c(F,T, T, T)
)

add_heatmap_track(Ath_Chromatin,Ath_Chromatin$group,fill=col_fun)

add_track(
    A,
    panel_fun = function(A) {
        x = A$start
        y = A$value
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "red",alpha = 0.5,lwd = 1))})

add_track(
    B,
    panel_fun = function(B) {
        x = B$start
        y = B$value
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "blue",alpha = 0.5,lwd = 1))},track = current_track())

add_rect_track(DLR, h1 = DLR$score, h2 = 0,
               gp = gpar(col = ifelse(DLR[[4]] > 0, "#EE0000FF", "#3B4992FF")))
add_rect_track(ICF, h1 = ICF$score, h2 = 0,
               gp = gpar(col = ifelse(ICF[[4]] > 0, "#EE0000FF", "#3B4992FF")))


dev.off()

## ----
Chr_df = data.frame(
    seqnames=c("Chr4"),
    start=c(1),
    end=c(6000000)
)

pdf(file = "UP_gtrellis_knob_ICF_DLR.pdf",height = 6,width = 10)
#range(UPTE[chr %in% c("Chr3","Chr4")]$value)
gtrellis_layout(
    data = Chr_df,
    gap = 0,
    track_ylim = c(c(0,1),c(-0.5,2.2), c(-0.05, 0.05), c(-0.2, 0.2)),
    asist_ticks = T,
    add_name_track = T,
    name_fontsize = 12,
    #legend = lgd,
    name_track_fill = "white",
    track_ylab = c("","TE", "DLR", "ICF"),
    axis_label_fontsize = 8,
    lab_fontsize = 12,
    remove_chr_prefix = F,
    n_track = 4,
    track_height = c(0.2,1,1, 1),
    track_axis = c(F,T, T, T)
)

add_heatmap_track(Ath_Chromatin,Ath_Chromatin$group,fill=col_fun)

add_track(
    A,
    panel_fun = function(A) {
        x = A$start
        y = A$value
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "red",alpha = 0.5,lwd = 1))})

add_track(
    B,
    panel_fun = function(B) {
        x = B$start
        y = B$value
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "blue",alpha = 0.5,lwd = 1))},track = current_track())

add_rect_track(DLR, h1 = DLR$score, h2 = 0,
               gp = gpar(col = ifelse(DLR[[4]] > 0, "#EE0000FF", "#3B4992FF")))
add_rect_track(ICF, h1 = ICF$score, h2 = 0,
               gp = gpar(col = ifelse(ICF[[4]] > 0, "#EE0000FF", "#3B4992FF")))


dev.off()