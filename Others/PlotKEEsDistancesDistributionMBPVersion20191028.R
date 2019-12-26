#!/usr/bin/env Rscript
## Import some R packages ----
library(agricolae)
library(pacman)
suppressWarnings(suppressPackageStartupMessages(
    p_load(
        agricolae,
        data.table,
        stringr,
        dplyr,
        seriation,
        dendextend,
        dendsort,
        knitr,
        ComplexHeatmap,
        circlize,
        ggpubr,
        ggstatsplot,
        ggedit,
        ggsci,
        VIM
    )
))
source("/Users/sunlinhua/Downloads/theme_linhua.R")

theme_set(theme_linhua(base_size = 20,legend = "none"))
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
             aspect.ratio = 1/0.618
)

## Import FISH distance distribution table ----

## Pre from Raw data from Dr. Jing ----

K3K4_Control<-fread("~/Downloads/Control_K3K4.txt")
K3K4_Heat<-fread("~/Downloads/Heat_K3K4.txt")
K3K4_Recovery<-fread("~/Downloads/Recovery_K3K4.txt")
K7K8_Control<-fread("~/Downloads/Control_K7K8.txt")
K7K8_Heat<-fread("~/Downloads/Heat_K7K8.txt")
K7K8_Recovery<-fread("~/Downloads/Recovery_K7K8.txt")


K3K4_Control_Vec<-melt(K3K4_Control)[!is.na(melt(K3K4_Control)$value)]$value
K3K4_Heat_Vec<-melt(K3K4_Heat)[!is.na(melt(K3K4_Heat)$value)]$value
K3K4_Recovery_Vec<-melt(K3K4_Recovery)[!is.na(melt(K3K4_Recovery)$value)]$value
K7K8_Control_Vec<-melt(K7K8_Control)[!is.na(melt(K7K8_Control)$value)]$value
K7K8_Heat_Vec<-melt(K7K8_Heat)[!is.na(melt(K7K8_Heat)$value)]$value
K7K8_Recovery_Vec<-melt(K7K8_Recovery)[!is.na(melt(K7K8_Recovery)$value)]$value

rbind(data.table(KEE="K3K4", variable="Control",value=K3K4_Control_Vec),
      data.table(KEE="K3K4", variable="Heat",value=K3K4_Heat_Vec),
      data.table(KEE="K3K4", variable="Recovery",value=K3K4_Recovery_Vec),
      data.table(KEE="K7K8", variable="Control",value=K7K8_Control_Vec),
      data.table(KEE="K7K8", variable="Heat",value=K7K8_Heat_Vec),
      data.table(KEE="K7K8", variable="Recovery",value=K7K8_Recovery_Vec)) -> kee



## KEE7/KEE8 ----

mtest<-kee[KEE=="K7K8"]
UU<-mtest[,boxplot.stats(value)$stats[4],by="variable"]
res<-kruskal(mtest$value,mtest$variable, group=TRUE, p.adj="bonferroni")$group
RES<-as.data.table(res,"variable")
RES<-RES[UU, on="variable"]
print(RES)

mtest[!is.na(mtest$value)] %>% group_by(variable) %>% summarise(n()) %>% as.data.table() -> N_Summary
colnames(N_Summary)<-c("Sample","Number")
kable(N_Summary)

pdf(file = "KEE7And8DistanceDistributionCompare.pdf",height = 8,width = 6)
ggplot(mtest[value<20], aes(x = variable, y = value)) +
    stat_boxplot(geom = 'errorbar', width = 0.3,color = c(pal_aaas()(9)[3], pal_aaas()(9)[2], pal_aaas()(9)[1])) +
    scale_fill_aaas(alpha = 1) +
    geom_jitter(size = 1, color = "lightgrey") +
    geom_boxplot(
        position = position_identity(),
        outlier.shape = NA,
        outlier.colour = "black",
        size = 0.5,
        #outlier.shape = 1,
        width = 0.4,
        outlier.size = NA,
        notch = T,fill=NA,
        notchwidth = 0,color=c(pal_aaas()(9)[3],pal_aaas()(9)[2],pal_aaas()(9)[1])
    ) +
   
    geom_text(
        data = RES,
        aes(x = variable, y = V1, label = groups),
        colour = "black",
        size = 6,
        vjust = 0,
        hjust = 0.5
    ) +
    labs(x = "", y = "Distance (μm)") -> P
    P+scale_y_continuous(limits = c(0, 10.0))
dev.off()

pdf(file = "KEE7And8DistanceDistributionCompareLimit6.pdf",
    height = 8,
    width = 6)
P+scale_y_continuous(limits = c(0, 6))
dev.off()

## KEE3/KEE4----

mtest<-kee[KEE=="K3K4"]
UU<-mtest[,boxplot.stats(value)$stats[4],by="variable"]
res<-kruskal(mtest$value,mtest$variable, group=TRUE, p.adj="bonferroni")$group
RES<-as.data.table(res,"variable")
RES<-RES[UU, on="variable"]
print(RES)

mtest[!is.na(mtest$value)] %>% group_by(variable) %>% summarise(n()) %>% as.data.table() -> N_Summary
colnames(N_Summary)<-c("Sample","Number")
kable(N_Summary)

pdf(file = "KEE3And4DistanceDistributionCompare.pdf",
    height = 8,
    width = 6)
ggplot(mtest[value < 20], aes(x = variable, y = value)) +
    stat_boxplot(geom = 'errorbar', width = 0.3,color = c(pal_aaas()(9)[3], pal_aaas()(9)[2], pal_aaas()(9)[1])) +
    scale_fill_aaas(alpha = 1) +
    geom_jitter(size = 1, color = "lightgrey") +
    geom_boxplot(
        position = position_identity(),
        outlier.shape = NA,
        outlier.colour = "black",
        size = 0.5,
        #outlier.shape = 1,
        width = 0.4,
        outlier.size = NA,
        notch = T,
        fill = NA,
        notchwidth = 0,
        color = c(pal_aaas()(9)[3], pal_aaas()(9)[2], pal_aaas()(9)[1])
    ) +
    geom_text(
        data = RES,
        aes(x = variable, y = V1, label = groups),
        colour = "black",
        size = 6,
        vjust = 0,
        hjust = 0.5
    ) +labs(x = "", y = "Distance (μm)") -> P
    P+scale_y_continuous(limits = c(0, 10.0))
dev.off()

pdf(file = "KEE3And4DistanceDistributionCompareLimit6.pdf",
    height = 8,
    width = 6)
    P+scale_y_continuous(limits = c(0, 6))
dev.off()