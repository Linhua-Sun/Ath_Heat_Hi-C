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

## Import sim tables
## Ath_Super<-fread("/data1/linhua/QIANLAB/PROJECT/ddccm/DMR_test/RUN/allc_2_meth/Ath_Superfamily_Sim.tsv")
## Ath_Super$count<-NULL

## Summary, count, and score----

UpTePlusTotal %>%
    group_by(group, Transposon_Family) %>%
    summarise(count = n()) %>%
    group_by(group) %>%
    mutate(freq = count / sum(count)) %>%
    ungroup() -> TeFamCountFreq

inner_join(
    reshape2::dcast(TeFamCountFreq, Transposon_Family ~ group, value.var = "count"),
    reshape2::dcast(TeFamCountFreq, Transposon_Family ~ group, value.var = "freq"),
    by = "Transposon_Family",
    suffix = c(".count", ".freq")
) %>% as.data.table() %$% .[, "Enrichment_Score" := Up_TE.freq / Total.freq][]  -> TeFamCountFreqEnrichScore

TeFamCountFreqEnrichScore[, sum_total := sum(Total.count)][, sum_up := sum(Up_TE.count, na.rm = T)][]

#TeFamCountFreqEnrichScore[!is.na(Enrichment_Score)] [order(Up_TE.count, decreasing = T)] %>% fwrite(file = "Heat_Activated_TE_Family_Enrichment_Analysis_20191124_TotalFilter.csv")
TeFamCountFreqEnrichScore[!is.na(Enrichment_Score)] [order(Up_TE.count, decreasing = T)] %>% fwrite(file = "Heat_Activated_TE_Family_Enrichment_Analysis_20191124.csv")
## Ploting ----
#pdf(
#    "Heat_Activated_TE_Family_Enrichment_Analysis_20191124_TotalFilter.pdf",
#    height = 6,
#    width = 12
#)
pdf(
    "Heat_Activated_TE_Family_Enrichment_Analysis_20191124.pdf",
    height = 6,
    width = 12
)
TeFamCountFreqEnrichScore[!is.na(Enrichment_Score)][Up_TE.count > 2][order(Enrichment_Score, decreasing = T)]  %>%
    ggbarplot(
        x = "Transposon_Family",
        y = "Enrichment_Score",
        order = BB[order(Enrichment_Score, decreasing = T)]$TE_Super_Family,
        width = 0.7,
        label = c("Up_TE.count"),
        lab.vjust = 1.5,
        lab.col = "white",
        fill = "black",
        alpha = 0.8,
        size = 0.5
    ) %>%
    ggpar(
        legend = "none",
        rotate = 0,
        xlab = "",
        ylab = "Enrichment score",
        font.x = c(15),
        x.text.angle = 45,
        font.y = c(15),
        font.tickslab = c(15),
        font.legend = c(15)
    ) +
    border() +
    geom_hline(yintercept = 1,
               lty = "dashed",
               color = "red") +
    theme(aspect.ratio = 1 / 3)
dev.off()