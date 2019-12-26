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
library(ggsci)

setwd("/data1/linhua/QIANLAB/PROJECT/hic/HOMER")
## Import data ----
rm(list=ls())

deltaPC1<-as.data.table(import("ByArmsHeatMinusControl_Res20K_Win60K_PC1.bedGraph"))
Control_PC1<-as.data.table(import("Control_PC1_R20KW60K_ByArms.PC1.bedGraph"))
Heat_PC1<-as.data.table(import("Heat_PC1_R20KW60K_ByArms.PC1.bedGraph"))


deltaPC1$score<-deltaPC1$score*(-1)
Control_PC1$score<-Control_PC1$score*(-1)
Heat_PC1$score<-Heat_PC1$score*(-1)

TRACK_LIMS<-c(c(0,1),c(-0.5,2.2),c(-3, 3), c(-3, 3), c(-1.8, 1.8))
Prefix<-"UpdateReByArmsHeatUpTE_"
## Set regions ----
Genome_df = data.frame(
    seqnames=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),
    start=c(1, 1, 1, 1, 1),
    end=c(30427671, 19698289, 23459830, 18585056, 26975502)
)

Chr12_df = data.frame(
    seqnames=c("Chr1", "Chr2"),
    start=c(1, 1),
    end=c(30427671, 19698289)
)

Chr34_df = data.frame(
    seqnames=c("Chr3", "Chr4"),
    start=c(1, 1),
    end=c(23459830, 18585056)
)


KNOB_df = data.frame(
    seqnames=c("Chr4"),
    start=c(1),
    end=c(6000000)
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

## All Chromosomes gtrellis plot  ----

pdf(file = paste(Prefix,"gtrellis_Genome_PC1.pdf",sep = ""),height = 6,width = 10)
gtrellis_layout(
    data = Genome_df,
    gap = 0,
    track_ylim = TRACK_LIMS,
    asist_ticks = T,
    add_name_track = T,
    name_fontsize = 12,
    name_track_fill = "white",
    track_ylab = c("","TE", "ControlPC1","HeatPC1","deltaPC1"),
    axis_label_fontsize = 8,
    lab_fontsize = 12,
    remove_chr_prefix = F,
    n_track = 5,
    track_height = c(0.2,1,1,1, 1),
    track_axis = c(F,T, T,T, T)
)

add_heatmap_track(Ath_Chromatin,Ath_Chromatin$group,fill=col_fun)

add_track(
    A,
    panel_fun = function(A) {
        x = A$start
        y = A$value
        grid.points(x,y,pch = 21,size = unit(2, "mm"),gp = gpar(col = "black",fill = "red",alpha = 0.5,lwd = 1))})

add_track(
    B,
    panel_fun = function(B) {
        x = B$start
        y = B$value
        grid.points(x,y,pch = 21,size = unit(2, "mm"),gp = gpar(col = "black",fill = "blue",alpha = 0.5,lwd = 1))},track = current_track())

add_rect_track(Control_PC1, h1 = Control_PC1$score, h2 = 0,
               gp = gpar(col = ifelse(Control_PC1$score > 0, "#EE0000FF", "#3B4992FF")))
add_rect_track(Heat_PC1, h1 = Heat_PC1$score, h2 = 0,
               gp = gpar(col = ifelse(Heat_PC1$score > 0, "#EE0000FF", "#3B4992FF")))
add_rect_track(deltaPC1, h1 = deltaPC1$score, h2 = 0,
               gp = gpar(col = ifelse(deltaPC1$score > 0, "#EE0000FF", "#3B4992FF")))
dev.off()

## Chromosomes3/4 gtrellis plot ----
pdf(file = paste(Prefix,"gtrellis_Chr34_PC1.pdf",sep = ""),height = 6,width = 10)
#range(UPTE[chr %in% c("Chr3","Chr4")]$value)
gtrellis_layout(
    data = Chr34_df,
    gap = 0,
    track_ylim = TRACK_LIMS,
    asist_ticks = T,
    add_name_track = T,
    name_fontsize = 12,
    #legend = lgd,
    name_track_fill = "white",
    track_ylab = c("","TE", "ControlPC1","HeatPC1","deltaPC1"),
    axis_label_fontsize = 8,
    lab_fontsize = 12,
    remove_chr_prefix = F,
    n_track = 5,
    track_height = c(0.2,1,1,1, 1),
    track_axis = c(F,T, T,T, T)
)

add_heatmap_track(Ath_Chromatin,Ath_Chromatin$group,fill=col_fun)

add_track(
    A,
    panel_fun = function(A) {
        x = A$start
        y = A$value
        grid.points(x,y,pch = 21,size = unit(2, "mm"),gp = gpar(col = "black",fill = "red",alpha = 0.5,lwd = 1))})

add_track(
    B,
    panel_fun = function(B) {
        x = B$start
        y = B$value
        grid.points(x,y,pch = 21,size = unit(2, "mm"),gp = gpar(col = "black",fill = "blue",alpha = 0.5,lwd = 1))},track = current_track())

add_rect_track(Control_PC1, h1 = Control_PC1$score, h2 = 0,
               gp = gpar(col = ifelse(Control_PC1$score > 0, "#EE0000FF", "#3B4992FF")))
add_rect_track(Heat_PC1, h1 = Heat_PC1$score, h2 = 0,
               gp = gpar(col = ifelse(Heat_PC1$score > 0, "#EE0000FF", "#3B4992FF")))
add_rect_track(deltaPC1, h1 = deltaPC1$score, h2 = 0,
               gp = gpar(col = ifelse(deltaPC1$score > 0, "#EE0000FF", "#3B4992FF")))
dev.off()
## ----
pdf(file = paste(Prefix,"gtrellis_KNOB_PC1.pdf",sep = ""),height = 6,width = 10)
#range(UPTE[chr %in% c("Chr3","Chr4")]$value)
gtrellis_layout(
    data = KNOB_df,
    gap = 0,
    track_ylim = TRACK_LIMS,
    asist_ticks = T,
    add_name_track = T,
    name_fontsize = 12,
    #legend = lgd,
    name_track_fill = "white",
    track_ylab = c("","TE", "ControlPC1","HeatPC1","deltaPC1"),
    axis_label_fontsize = 8,
    lab_fontsize = 12,
    remove_chr_prefix = F,
    n_track = 5,
    track_height = c(0.2,1,1,1, 1),
    track_axis = c(F,T, T,T, T)
)

add_heatmap_track(Ath_Chromatin,Ath_Chromatin$group,fill=col_fun)

add_track(
    A,
    panel_fun = function(A) {
        x = A$start
        y = A$value
        grid.points(x,y,pch = 21,size = unit(2, "mm"),gp = gpar(col = "black",fill = "red",alpha = 0.5,lwd = 1))})

add_track(
    B,
    panel_fun = function(B) {
        x = B$start
        y = B$value
        grid.points(x,y,pch = 21,size = unit(2, "mm"),gp = gpar(col = "black",fill = "blue",alpha = 0.5,lwd = 1))},track = current_track())

add_rect_track(Control_PC1, h1 = Control_PC1$score, h2 = 0,
               gp = gpar(col = ifelse(Control_PC1$score > 0, "#EE0000FF", "#3B4992FF")))
add_rect_track(Heat_PC1, h1 = Heat_PC1$score, h2 = 0,
               gp = gpar(col = ifelse(Heat_PC1$score > 0, "#EE0000FF", "#3B4992FF")))
add_rect_track(deltaPC1, h1 = deltaPC1$score, h2 = 0,
               gp = gpar(col = ifelse(deltaPC1$score > 0, "#EE0000FF", "#3B4992FF")))
dev.off()
