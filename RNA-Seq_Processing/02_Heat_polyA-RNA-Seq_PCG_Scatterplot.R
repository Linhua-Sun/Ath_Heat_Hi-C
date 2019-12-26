## Library R packages ----
library(pacman)
suppressWarnings(suppressPackageStartupMessages(
    p_load(
        ComplexHeatmap,
        DESeq2,
        crayon,
        data.table,
        dendsort,
        dplyr,
        edgeR,
        ggplot2,
        ggpubr,
        ggsci,
        gmodels,
        hashmap,
        readr,
        rtracklayer,
        scatterplot3d,
        stringr,
        tidyverse,
        tximport,
        ggrepel,
        viridis
    )
))

source("/data1/linhua/QIANLAB/R_Functions/theme_linhua.R")

theme_set(theme_linhua(base_size = 18,base_family = "",legend = "right"))
theme_update(panel.grid.minor.x = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.grid.major.x = element_blank(),
             panel.grid.major.y = element_blank(),
             #axis.text.x = element_blank(),
             #axis.title.x = element_blank(),
             panel.spacing.x = unit(0.5, "lines"),
             panel.spacing.y = unit(0.5, "lines"),
             strip.background = element_blank(),
             panel.border = element_rect(colour = "black"),
             aspect.ratio = 1
)

## DEG (GENE) Pairwise compare scatterplot ----

GENE<-fread("Updated_Heat_Stress_polyA-RNA-Seq_DEG_Ann_Results.csv")

## Choose selected genes highlight in scatterplot ----
HSP_ID<-c("AT1G16030","AT1G74310","AT5G12030","AT3G13860","AT4G25200","AT3G46230","AT5G52640","AT3G07770","AT4G27670","AT5G12020","AT3G23990","AT3G12580","AT5G59720")
HSP_Name<-c("Hsp70b","HSP101","HSP17.6A","HSP60-3A","HSP23.6-MITO","HSP17.4","HSP90.1","Hsp89.1","HSP21","HSP17.6II","HSP60","HSP70","HSP18.2")
HSP <-data.table(ID=HSP_ID,Name=toupper(HSP_Name))

GENE[ID=="AT3G30720"] # QQS
GENE[ID=="AT2G17690"] # SDC
GENE[ID=="AT2G36490"]$symbol<-"ROS1"
GENE[ID=="AT1G06760"]$symbol<-"H1.1"
GENE[ID=="AT2G30620"]$symbol<-"H1.2"
GENE[ID=="AT2G18050"]$symbol<-"H1.3"

DOWN_Marker<-GENE[ID %in% c("AT2G36490","AT1G06760","AT2G30620","AT2G18050","AT3G30720","AT2G17690")][,c("ID","symbol")]
setnames(DOWN_Marker,"symbol","Name")

SE <-
    HSP[Name %in% toupper(
        c(
            "HSP90.1",
            "HSP17.6II",
            "HSP18.2",
            "HSP21",
            "HSP70",
            "HSP101",
            "HSP17.6A",
            "HSP17.4",
            "HSP23.6−MITO"
        )
    )]

SE<-rbind(SE,DOWN_Marker)

## Use CPM to plot scatterplot ----
P <-ggplot(GENE, aes(x = Control.mean.cpm.log2, y = Heat.mean.com.log2)) +
    #geom_point(
    #    data = GENE[Heat_vs_Control == "NoChange"],
    #    aes(x = Control.mean.cpm.log2, y = Heat.mean.com.log2),
    #    color = "black",
    #    alpha = 0.5,
    #    size = 1
    #) +
    scale_x_continuous(limits = c(0, 15), name = "Control (log2(cpm))") +
    scale_y_continuous(limits = c(0, 15), name = "Heat (log2(cpm))")

L1 <- P + stat_density_2d(aes(fill = ..level..), geom = "polygon") +
    scale_fill_viridis(option = "D")

L2 <- geom_point(data = GENE[Heat_vs_Control=="UP"],aes(x = Control.mean.cpm.log2, y = Heat.mean.com.log2),color = "#BE2442",alpha = 0.4,size = 2)

L3 <- geom_point(data = GENE[Heat_vs_Control=="DOWN"],aes(x = Control.mean.cpm.log2, y = Heat.mean.com.log2),color = "#324743",alpha = 0.4,size = 2)

L4_highlight<-geom_point(data = GENE[ID %in% SE$ID],aes(x = Control.mean.cpm.log2, y = Heat.mean.com.log2),color = "black",alpha = 1,size = 3,fill = 'green',shape = 24)

L4_label<-geom_text_repel(data = GENE[ID %in% SE$ID],aes(x = Control.mean.cpm.log2, y = Heat.mean.com.log2,label=GENE[ID %in% SE$ID]$symbol),color="blue",alpha=1,size=5)

L1 +L2 + L3 +L4_highlight +L4_label+
    geom_abline(slope = 1, intercept = log2(2),lty="solid",color = "#FF8000",size = 1)+
    geom_abline(slope = 1,intercept = -log2(2),lty="solid",color = "#FF8000",size = 1) -> GENES_Scatterplot_CPM
#geom_abline(intercept = 0,color = "black",size = 0.5)

## Use FPKM to plot scatterplot ----

pdf(file = "Heat_GENE_Pairwise_Compare_SE_HSPs_RPKM_Detailed.pdf",height = 6,width = 6)

P <-ggplot(GENE, aes(x = Control.mean.rpkm.log2, y = Heat.mean.rpkm.log2)) +
    scale_x_continuous(limits = c(0, 15), name = "Control (log2(FPKM))") +
    scale_y_continuous(limits = c(0, 15), name = "Heat (log2(FPKM))")

L1 <- P + stat_density_2d(aes(fill = ..level..), geom = "polygon") +
    scale_fill_viridis(option = "D")

L2 <- geom_point(data = GENE[Heat_vs_Control=="UP"],aes(x = Control.mean.rpkm.log2, y = Heat.mean.rpkm.log2),color = "#BE2442",alpha = 0.4,size = 2)

L3 <- geom_point(data = GENE[Heat_vs_Control=="DOWN"],aes(x = Control.mean.rpkm.log2, y = Heat.mean.rpkm.log2),color = "#324743",alpha = 0.4,size = 2)

L4_highlight<-geom_point(data = GENE[ID %in% SE$ID],aes(x = Control.mean.rpkm.log2, y = Heat.mean.rpkm.log2),color = "black",alpha = 1,size = 3,fill = 'green',shape = 24)

L4_label<-geom_text_repel(data = GENE[ID %in% SE$ID],aes(x = Control.mean.rpkm.log2, y = Heat.mean.rpkm.log2,label=GENE[ID %in% SE$ID]$symbol),color="blue",alpha=1,size=5)

L1 +L2 + L3 +L4_highlight +L4_label+
    geom_abline(slope = 1, intercept = log2(2),lty="solid",color = "#FF8000",size = 1)+
    geom_abline(slope = 1,intercept = -log2(2),lty="solid",color = "#FF8000",size = 1) 

dev.off()

GENE[ID %in% SE$ID][,c("ID","Heat_vs_Control","Recovery_vs_Heat","symbol","Control.mean.rpkm.log2","Heat.mean.rpkm.log2","Recovery.mean.rpkm.log2","Control.mean.rpkm","Heat.mean.rpkm","Recovery.mean.rpkm","Control.mean.cpm","Heat.mean.cpm","Recovery.mean.cpm")][order(symbol)] %>% kable(format = "rst",digits = 0) %>%  write(file = "Genes_Scatterplot_Highlight_IDs.txt")