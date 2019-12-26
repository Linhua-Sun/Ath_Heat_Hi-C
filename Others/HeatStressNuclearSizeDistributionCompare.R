#!/usr/bin/env Rscript
## Import some R packages ----

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

CC<-fread("~/Desktop/22_nuclear_size.csv")[,"Group":="Control"]
TT<-fread("~/Desktop/37_nuclear_size.csv")[,"Group":="Heat"]

library(agricolae)

mtest<-rbind(CC,TT)
colnames(mtest)<-c("ID","value","variable")
##
##UU<-mtest[,boxplot.stats(value)$stats[4],by="variable"]
##res<-kruskal(mtest$value,mtest$variable, group=TRUE, p.adj="bonferroni")$group
##RES<-as.data.table(res,"variable")
##RES<-RES[UU, on="variable"]
##print(RES)
##
##mtest[!is.na(mtest$value)] %>% group_by(variable) %>% summarise(n()) %>% as.data.table() -> N_Summary
##colnames(N_Summary)<-c("Sample","Number")
##kable(N_Summary)

pdf(file = "HeatStressNuclearSizedistributionCompare.pdf",height = 10,width = 6)

ggplot(mtest, aes(x = variable, y = value,fill= variable)) +
    geom_boxplot(
        position = position_identity(),
        outlier.shape = NA,
        outlier.colour = "black",
        size = 0.5,
        #outlier.shape = 1,
        width = 0.4,
        outlier.size = NA,
        notch = T,
        notchwidth = 0,fill=c(pal_aaas()(9)[3],pal_aaas()(9)[2])
    ) +
    stat_boxplot(geom = 'errorbar', width = 0.3) +
    geom_jitter(size = 1, color = "darkgrey") +
    #scale_fill_aaas(alpha = 1) +
    #geom_text(
    #    data = RES,
    #    aes(x = variable, y = V1, label = groups),
    #    colour = "red",
    #    size = 6,
    #    vjust = 0,
    #    hjust = 0.5
    #) +
    labs(x = "", y = "Diameter (μm)")
dev.off()

compare_means(formula = value ~ variable,data = mtest,method ="wilcox.test" )
