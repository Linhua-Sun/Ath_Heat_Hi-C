TT<-fread("Updated_Heat_Stress_polyA-RNA-Seq_DEG_TEs_Ann_Results.csv")
UP<-TT[Heat_vs_Control=="UP"]

## ----

library(gtrellis)

#UP$log2_FC<-log2(UP$Heat.mean.com.log2 / (UP$Control.mean.cpm.log2+0.001))

#UP[log2_FC>8]$log2_FC<-8

UP$log10<-log10(UP$Heat.mean.cpm - UP$Control.mean.cpm)

GRUP<-GRanges(UP)

Ath_chr = fread("/data1/linhua/software/BIN/chrom.sizes")[1:5,]
Chr_df = data.frame(c("Chr1","Chr2","Chr3", "Chr4", "Chr5"), c(1,1,1,1,1),c(30427671,19698289, 23459830, 18585056, 26975502))
centromere=data.table(chr=paste0("Chr",1:5),start=c(13000000,2000000,11500000,3000000,10000000),end=c(17000000,5000000,15500000,5000000,14000000),val=1)

Ath_Chromatin<-fread("/data1/linhua/QIANLAB/PROJECT/hic/HOMER/BEDG/Ath_cytoband_Centromere-200K_gieStain.csv")


Ath_Chromatin$group<-0
Ath_Chromatin[gieStain=="acen"]$group<-1
Ath_Chromatin[gieStain=="gneg"]$group<-2
Ath_Chromatin[gieStain=="gpos100"]$group<-3
Ath_Chromatin[gieStain=="gpos95"]$group<-4
Ath_Chromatin[gieStain=="stalk"]$group<-5


#TE<-fread("/data1/linhua/QIANLAB/PROJECT/ATAC-Seq/ATAC-seq/ATAC_Out/Reads_Counts/ME/sort_TAIR10_TE.bed",showProgress = T)
#TE_density<-circlize::genomicDensity(region = TE,window.size = 1e6)
#TE_density$chr<-str_replace_all(TE_density$chr,"chr","Chr")
col_fun = circlize::colorRamp2(seq(1,5),c("red","white","black","grey","red"))



Chr_df = data.frame(c("Chr1","Chr2","Chr3", "Chr4", "Chr5"), c(1,1,1,1,1),c(30427671,19698289, 23459830, 18585056, 26975502))
Chr_df1 <- data.frame(c("Chr1"), c(1),c(30427671))
Chr_df12 <- data.frame(c("Chr1","Chr2"), c(1,1),c(30427671,19698289))
Chr_df3 <- data.frame(c("Chr3"), c(1),c(23459830))
Chr_df4 <- data.frame(c("Chr4"), c(1),c(18585056))

Chr_KNOB <- data.frame(c("Chr4"), c(1),c(6000000))

Chr_df34 <- data.frame(c("Chr3","Chr4"), c(1,1),c(23459830,18585056))




GG<-GRUP[GRUP$Transposon_Super_Family=="LTR/Gypsy"]
CC<-GRUP[GRUP$Transposon_Super_Family=="LTR/Copia"]
OT<-GRUP[!(GRUP$Transposon_Super_Family %in% c("LTR/Copia","LTR/Gypsy")) & !(GRUP$Class %in% "DNAtransposon")]
DNA<-GRUP[GRUP$Class=="DNAtransposon"]

pdf("TE_distr_in_Chr_df12.pdf",height = 6,width = 12)

gtrellis_layout(data=Chr_df12,
                n_track=3,
                #track_ylab=c("",""),
                asist_ticks=F,
                title='TE distribution in chromosome 1 and 2',
                track_ylim=c(c(-0.5,2.5),c(-0.5,2.5),c(0,1)),
                track_height=c(2,2,0.3),
                track_axis=c(T,T,F))
add_track(
    OT,
    panel_fun = function(OT) {
        x = start(OT)
        y = OT$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "grey",alpha = 0.7,lwd = 1))})


add_track(
    CC,
    panel_fun = function(CC) {
        x = start(CC)
        y = CC$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "blue",alpha = 0.7,lwd = 1))},track = current_track())

add_track(
    GG,
    panel_fun = function(GG) {
        x = start(GG)
        y = GG$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "red",alpha = 0.5,lwd = 1))},track = current_track())


add_track(
    DNA,
    panel_fun = function(DNA) {
        x = start(DNA)
        y = DNA$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "green",alpha = 0.7,lwd = 1))})

add_heatmap_track(Ath_Chromatin,Ath_Chromatin$group,fill=col_fun)

dev.off()

## Genome-Wide TE distribution ----

TE_Annotation <- fread("/data1/linhua/QIANLAB/PROJECT/te/TE_Annotation-2018-11-17.csv")
table(TE_Annotation$Transposon_Super_Family)
COPIA<-GRanges(TE_Annotation[Transposon_Super_Family=="LTR/Copia"])

Copia_density<-circlize::genomicDensity(region = TE_Annotation[Transposon_Super_Family=="LTR/Copia"][,c("seqnames","start","end")],window.size = 1e5,overlap = F)

Gypsy_Density<-circlize::genomicDensity(region = TE_Annotation[Transposon_Super_Family=="LTR/Gypsy"][,c("seqnames","start","end")],window.size = 1e5,overlap = F)

DNA_TE_Density<-circlize::genomicDensity(region = TE_Annotation[Class=="DNAtransposon"][,c("seqnames","start","end")],window.size = 1e5,overlap = F)

## TE_distr_in_Chr12 ----
pdf("TE_distr_in_Chr12.pdf",height = 6,width = 12)

gtrellis_layout(data=Chr_df12,
                n_track=4,
                #track_ylab=c("",""),
                asist_ticks=F,
                title='TE distribution in Chr 1 2',
                track_ylim=c(c(0,1),range(GRUP$log10),range(GRUP$log10),c(0,1)),
                track_height=c(0.5,2,2,0.3),
                track_axis=c(T,T,T,F))
add_lines_track(gr = DNA_TE_Density,
                value = DNA_TE_Density$pct/max(DNA_TE_Density$pct),
                gp = gpar(col = "green"))

add_lines_track(
    gr = Copia_density,
    value = Copia_density$pct,
    gp = gpar(col = "blue"),
    track = current_track()
)

add_lines_track(
    gr = Gypsy_Density,
    value = Gypsy_Density$pct,
    gp = gpar(col = "red"),
    track = current_track()
)



add_track(
    OT,
    panel_fun = function(OT) {
        x = start(OT)
        y = OT$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "grey",alpha = 0.7,lwd = 1))})


add_track(
    CC,
    panel_fun = function(CC) {
        x = start(CC)
        y = CC$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "blue",alpha = 0.7,lwd = 1))},track = current_track())

add_track(
    GG,
    panel_fun = function(GG) {
        x = start(GG)
        y = GG$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "red",alpha = 0.5,lwd = 1))},track = current_track())


add_track(
    DNA,
    panel_fun = function(DNA) {
        x = start(DNA)
        y = DNA$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "green",alpha = 0.7,lwd = 1))})

add_heatmap_track(Ath_Chromatin,Ath_Chromatin$group,fill=col_fun)

dev.off()
## TE_distr_in_Genome.pdf ----
pdf("TE_distr_in_Genome.pdf",height = 6,width = 12)

gtrellis_layout(data=Chr_df,
                n_track=4,
                #track_ylab=c("",""),
                asist_ticks=F,
                title='TE distribution in Chr 1 2',
                track_ylim=c(c(0,1),range(GRUP$log10),range(GRUP$log10),c(0,1)),
                track_height=c(0.5,2,2,0.3),
                track_axis=c(T,T,T,F))
add_lines_track(gr = DNA_TE_Density,
                value = DNA_TE_Density$pct/max(DNA_TE_Density$pct),
                gp = gpar(col = "green"))

add_lines_track(
    gr = Copia_density,
    value = Copia_density$pct,
    gp = gpar(col = "blue"),
    track = current_track()
)

add_lines_track(
    gr = Gypsy_Density,
    value = Gypsy_Density$pct,
    gp = gpar(col = "red"),
    track = current_track()
)



add_track(
    OT,
    panel_fun = function(OT) {
        x = start(OT)
        y = OT$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "grey",alpha = 0.7,lwd = 1))})


add_track(
    CC,
    panel_fun = function(CC) {
        x = start(CC)
        y = CC$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "blue",alpha = 0.7,lwd = 1))},track = current_track())

add_track(
    GG,
    panel_fun = function(GG) {
        x = start(GG)
        y = GG$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "red",alpha = 0.5,lwd = 1))},track = current_track())


add_track(
    DNA,
    panel_fun = function(DNA) {
        x = start(DNA)
        y = DNA$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "green",alpha = 0.7,lwd = 1))})

add_heatmap_track(Ath_Chromatin,Ath_Chromatin$group,fill=col_fun)

dev.off()

## TE_distr_in_Chr34 ----
pdf("TE_distr_in_Chr34.pdf",height = 6,width = 12)

gtrellis_layout(data=Chr_df34,
                n_track=4,
                #track_ylab=c("",""),
                asist_ticks=F,
                title='TE distribution in Chr 1 2',
                track_ylim=c(c(0,1),range(GRUP$log10),range(GRUP$log10),c(0,1)),
                track_height=c(0.5,2,2,0.3),
                track_axis=c(T,T,T,F))
add_lines_track(gr = DNA_TE_Density,
                value = DNA_TE_Density$pct/max(DNA_TE_Density$pct),
                gp = gpar(col = "green"))

add_lines_track(
    gr = Copia_density,
    value = Copia_density$pct,
    gp = gpar(col = "blue"),
    track = current_track()
)

add_lines_track(
    gr = Gypsy_Density,
    value = Gypsy_Density$pct,
    gp = gpar(col = "red"),
    track = current_track()
)



add_track(
    OT,
    panel_fun = function(OT) {
        x = start(OT)
        y = OT$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "grey",alpha = 0.7,lwd = 1))})


add_track(
    CC,
    panel_fun = function(CC) {
        x = start(CC)
        y = CC$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "blue",alpha = 0.7,lwd = 1))},track = current_track())

add_track(
    GG,
    panel_fun = function(GG) {
        x = start(GG)
        y = GG$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "red",alpha = 0.5,lwd = 1))},track = current_track())


add_track(
    DNA,
    panel_fun = function(DNA) {
        x = start(DNA)
        y = DNA$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "green",alpha = 0.7,lwd = 1))})

add_heatmap_track(Ath_Chromatin,Ath_Chromatin$group,fill=col_fun)

dev.off()
## TE_distr_in_KNOB ----
pdf("TE_distr_in_KNOB.pdf",height = 6,width = 12)

gtrellis_layout(data=Chr_KNOB,
                n_track=4,
                #track_ylab=c("",""),
                asist_ticks=F,
                title='TE distribution in Chr 1 2',
                track_ylim=c(c(0,1),range(GRUP$log10),range(GRUP$log10),c(0,1)),
                track_height=c(0.5,2,2,0.3),
                track_axis=c(T,T,T,F))
add_lines_track(gr = DNA_TE_Density,
                value = DNA_TE_Density$pct/max(DNA_TE_Density$pct),
                gp = gpar(col = "green"))

add_lines_track(
    gr = Copia_density,
    value = Copia_density$pct,
    gp = gpar(col = "blue"),
    track = current_track()
)

add_lines_track(
    gr = Gypsy_Density,
    value = Gypsy_Density$pct,
    gp = gpar(col = "red"),
    track = current_track()
)



add_track(
    OT,
    panel_fun = function(OT) {
        x = start(OT)
        y = OT$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "grey",alpha = 0.7,lwd = 1))})


add_track(
    CC,
    panel_fun = function(CC) {
        x = start(CC)
        y = CC$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "blue",alpha = 0.7,lwd = 1))},track = current_track())

add_track(
    GG,
    panel_fun = function(GG) {
        x = start(GG)
        y = GG$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "red",alpha = 0.5,lwd = 1))},track = current_track())


add_track(
    DNA,
    panel_fun = function(DNA) {
        x = start(DNA)
        y = DNA$log10
        grid.points(x,y,pch = 21,size = unit(3, "mm"),gp = gpar(col = "black",fill = "green",alpha = 0.7,lwd = 1))})

add_heatmap_track(Ath_Chromatin,Ath_Chromatin$group,fill=col_fun)

dev.off()