## ----
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(gtrellis))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(ggsci))

## ----
Prefix<-"Heat_Activated_TE_SuperFamily_Enrichment"
Ath_Super<-fread("/data1/linhua/QIANLAB/PROJECT/ddccm/DMR_test/RUN/allc_2_meth/Ath_Superfamily_Sim.tsv")
Ath_Super$count<-NULL

## Import TE expression and TE annotation ----

TE_Annotation <-
    fread(
        "/data1/linhua/QIANLAB/PROJECT/te/TE_Annotation-2019-11-21-PlusH1AndMNasePlusLocation.csv"
    )

#TE_Annotation<-TE_Annotation[!is.na(score.GC)] ## _TotalFilter

TE_ExprAnn <-
    fread(
        "/data1/linhua/QIANLAB/PROJECT/te/HeatTEReanalysis2019-06-03/Updated_Heat_Stress_polyA-RNA-Seq_DEG_TEs_Ann_Results.csv"
    )

TE_ExprAnn[Heat_vs_Control == "UP"][, intersect(colnames(TE_Annotation), colnames(TE_ExprAnn)), with =
                                        F][, "group" := "Up_TE"] -> Up_TE_ANN
TE_Annotation[, intersect(colnames(TE_Annotation), colnames(TE_ExprAnn)), with =
                  F][, "group" := "Total"] -> Total_ANN

rbind(Up_TE_ANN, Total_ANN) -> UpTePlusTotal
UpTePlusTotal$group <- factor(UpTePlusTotal$group)

## Similifyed TE class enrichment analysis----
UpTePlusTotal %>% group_by(group,Class) %>% summarise(count=n()) %>% mutate(freq=count/sum(count)) %>% fwrite(paste(Prefix,"_TE_Class_Simiple.csv",sep = ""))
pdf(file = paste(Prefix,"_TE_Class_Simiple.pdf",sep = ""),width = 8,height = 8)
UpTePlusTotal %>% group_by(group,Class) %>% summarise(count=n()) %>% mutate(freq=count/sum(count)) %>%
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
    theme(aspect.ratio = 3,axis.line.x =element_blank(),axis.ticks.x = element_blank())
dev.off()


## ----
left_join(UpTePlusTotal, Ath_Super, by = "Transposon_Super_Family") %>%
    group_by(group, TE_Super_Family) %>%
    summarise(count = n()) %>%
    group_by(group) %>%
    mutate(freq = count / sum(count)) %>%
    ungroup() -> TeSuperFamCountFreq

inner_join(
    reshape2::dcast(TeSuperFamCountFreq, TE_Super_Family ~ group, value.var = "count"),
    reshape2::dcast(TeSuperFamCountFreq, TE_Super_Family ~ group, value.var = "freq"),
    by = "TE_Super_Family",
    suffix = c(".count", ".freq")
) %>% as.data.table() %$% .[TE_Super_Family != "Other"][, "Enrichment_Score" := Up_TE.freq / Total.freq][]  -> TeFamCountFreqEnrichScore

TeFamCountFreqEnrichScore[, sum_total := sum(Total.count)][, sum_up := sum(Up_TE.count)][]


## Finnal barplot ----
TeFamCountFreqEnrichScore %>% fwrite(paste(Prefix,"_TE_SuperFamily_Enrichment.csv",sep = ""))
pdf(file = paste(Prefix,"_TE_SuperFamily_Enrichment.pdf",sep = ""),
       width = 12,
       height = 6)
TeFamCountFreqEnrichScore  %>%
    ggbarplot(
        x = "TE_Super_Family",
        y = "Enrichment_Score",
        #fill = "condition",
        order = TeFamCountFreqEnrichScore[order(Enrichment_Score, decreasing = T)]$TE_Super_Family,
        width = 0.7,
        #fill = pal_aaas()(2)[2],
        label = c("Up_TE.count"),
        lab.vjust = 1.5,
        lab.col = "white",
        fill = "black",
        alpha = 0.8,
        size = 0.5
    ) %>%
    ggpar(
        legend = "none",
        xlab = "",
        ylab = "Enrichment score",
        font.x = c(15),
        x.text.angle = 45,
        font.y = c(15),
        font.tickslab = c(15),
        font.legend = c(15)
    ) +
    border() +
    #scale_x_discrete(expand = c(0, 0)) +
    #scale_y_continuous(expand = c(0, 0))
    geom_hline(yintercept = 1,
               lty = "dashed",
               color = "blue") +
    #theme(aspect.ratio = (1 + sqrt(5)) / 2) # golden ratio portrait
    theme(aspect.ratio = 0.618)
    #geom_text(data = temp,aes(x = TE_Super_Family, y = Enrichment_Score,label= format),color="red",size=10)
dev.off()


## barplot similar to GB 2016 ----
reshape2::dcast(TeSuperFamCountFreq, TE_Super_Family ~ group,value.var="freq") %>% as.data.table() -> LLL
LLL[TE_Super_Family!="Other"][,"log2":=log2(Up_TE/Total)][]  -> AA

AA$condition<-"*"
AA[log2>0]$condition<-"UP"
AA[log2<0]$condition<-"DOWN"

pdf(file = paste(Prefix,"_TE_SuperFamily_Enrichment_BidirectionalBars.pdf",sep = ""),width = 12,height = 6)
AA %>%
    ggbarplot(
        x = "TE_Super_Family",
        y = "log2",
        fill = "condition",
        width = 0.7,
        palette = "aaas",
        alpha = 0.8,
        size = 0.5
    ) %>%
    ggpar(
        legend = "none",
        xlab = "",
        ylab = "Enrichment",
        font.x = c(18),
        x.text.angle = 30,
        font.y = c(18),
        font.tickslab = c(15),
        font.legend = c(15)
    ) +
    #border() +
    #scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0),limits = c(-3,2)) + # https://stackoverflow.com/questions/13701347/force-the-origin-to-start-at-0
    geom_hline(yintercept = 0) +
    #theme(aspect.ratio = (1 + sqrt(5)) / 2) # golden ratio portrait
    theme(
        aspect.ratio = 0.8,
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()
    )
dev.off()

## Detailed TE superfamily enrichment analysis ----


pdf(file = paste(Prefix,"_Detailled_TE_Frac_2019.pdf",sep=""),width = 8,height = 8)

left_join(UpTePlusTotal,Ath_Super,by="Transposon_Super_Family") %>%
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