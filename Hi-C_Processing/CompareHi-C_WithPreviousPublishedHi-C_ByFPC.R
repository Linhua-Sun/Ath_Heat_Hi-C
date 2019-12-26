## Import R packages and Define R functions ----
library(pacman)
suppressWarnings(suppressPackageStartupMessages(
    p_load(
        HiCdatR,
        data.table,
        ggplot2,
        scales,
        ggpmisc,
        ggsci,
        ggpubr,
        sitools,
        extrafont,
        hashmap,
        grid,
        magrittr,
        dplyr
    )
))
suppressMessages(loadfonts())

source("/data1/linhua/QIANLAB/R_Functions/theme_linhua.R")

theme_set(theme_linhua(base_size = 20,base_family = "Arial",legend = "right",x.text.angle=0))

theme_update(#panel.grid.minor.x = element_blank(),
    #panel.grid.minor.y = element_blank(),
    #panel.grid.major.x = element_blank(),
    #panel.grid.major.y = element_blank(),
    #axis.text.x = element_blank(),
    #axis.title.x = element_blank(),
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(0.5, "lines"),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    aspect.ratio = 2 / (1 + sqrt(5)),
    #aspect.ratio = 1,
    plot.title = element_text(hjust=0, size=18),
    legend.title=element_blank(),
    axis.text.y = element_text(hjust = 1),
    axis.text.x = element_text(vjust=1)
)

base_breaks <- function(n = 5){
    function(x) {
        axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
    }
}

set_font <- function(p, family="sans", fontface=NULL, size=NULL, color=NULL) {
    library(grid)
    library(ggplot2)
    if (!is.null(size))
        size <- size  * .pt
    par <- list(fontfamily = family, fontface = fontface, fontsize = size, col = color)
    par <- par[!sapply(par, is.null)]
    gp <- do.call(gpar, par)
    g <- ggplotGrob(p)
    ng <- grid.ls(grid.force(g), print=FALSE)$name
    txt <- ng[which(grepl("text", ng))]
    
    for (i in seq_along(txt)) {
        g <- editGrob(grid.force(g), gPath(txt[i]),
                      grep = TRUE, gp = gp)
    }
    return(g)
}


## Import HiC Contact Matrix from HiC-Pro ----

setwd("/data1/linhua/QIANLAB/PROJECT/hic/FPC/BIN100K_Col_Matrix/BU")
pathToTutorial <- "/data1/linhua/QIANLAB/PROJECT/hic/FPC/BIN100K_Col_Matrix"
pathToScripts <- "/data1/linhua/QIANLAB/PROJECT/hic/MWSchmid_test/Rscripts/"
f.source.organism.specific.code(file.path(pathToScripts, "HiCdat-A-thaliana-TAIR10.R"))
definedGenomicRegionsArms <- read.table(file.path("/data1/linhua/QIANLAB/PROJECT/hic/MWSchmid_test/At_tutorial_files", "genomicRegionsArms_forPCA.txt"), sep = '\t', header = TRUE, quote = "", stringsAsFactors = FALSE)
binSize <- 100000

sampleList <- list(
    Col22_Qian_Mboi= c("Col22_total_sample_100000.matrix"),
    Col_2012_hindiii= c("Col_2012_hindiii_100000.matrix"),
    Col_Feng_hindiii= c("Col_Feng_hindiii_100000.matrix"),
    Col_Grob_hindiii= c("Col_Grob_hindiii_100000.matrix"),
    Col_Chang_dpnii=c("Col_Chang_dnpii_100000.matrix"),
    Col_Chang_dpnii_2018=c("Col_Chang_dpnii_2018_100000.matrix")
)

binMatList <- f.load.samples(pathToTutorial, sampleList, binSize, 50)

Col22_Qian_Mboi_PCA <- f.principle.component.analysis.and.features(dataMatrix = binMatList$Col22_Qian_Mboi,outfilePrefix = "Col22_Qian_Mboi_PCA",binSize = 100000,rDir = "./",regionTable = definedGenomicRegionsArms,simplifiedNames = list(),filterZero = TRUE,filterThreshold = 0.95,pValueThreshold = 0.05,userLimits = c(-1, 1))
Col_2012_hindiii_PCA <- f.principle.component.analysis.and.features(dataMatrix = binMatList$Col_2012_hindiii,outfilePrefix = "Col_2012_hindiii_PCA",binSize = 100000,rDir = "./",regionTable = definedGenomicRegionsArms,simplifiedNames = list(),filterZero = TRUE,filterThreshold = 0.95,pValueThreshold = 0.05,userLimits = c(-1, 1))
Col_Feng_hindiii_PCA <- f.principle.component.analysis.and.features(dataMatrix = binMatList$Col_Feng_hindiii,outfilePrefix = "Col_Feng_hindiii_PCA",binSize = 100000,rDir = "./",regionTable = definedGenomicRegionsArms,simplifiedNames = list(),filterZero = TRUE,filterThreshold = 0.95,pValueThreshold = 0.05,userLimits = c(-1, 1))
Col_Grob_hindiii_PCA <- f.principle.component.analysis.and.features(dataMatrix = binMatList$Col_Grob_hindiii,outfilePrefix = "Col_Grob_hindiii_PCA",binSize = 100000,rDir = "./",regionTable = definedGenomicRegionsArms,simplifiedNames = list(),filterZero = TRUE,filterThreshold = 0.95,pValueThreshold = 0.05,userLimits = c(-1, 1))

Col_Chang_dpnii_PCA <- f.principle.component.analysis.and.features(dataMatrix = binMatList$Col_Chang_dpnii,outfilePrefix = "Col_Chang_dpnii_PCA",binSize = 100000,rDir = "./",regionTable = definedGenomicRegionsArms,simplifiedNames = list(),filterZero = TRUE,filterThreshold = 0.95,pValueThreshold = 0.05,userLimits = c(-1, 1))
Col_Chang_dpnii_2018_PCA <- f.principle.component.analysis.and.features(dataMatrix = binMatList$Col_Chang_dpnii_2018,outfilePrefix = "Col_Chang_dpnii_2018_PCA",binSize = 100000,rDir = "./",regionTable = definedGenomicRegionsArms,simplifiedNames = list(),filterZero = TRUE,filterThreshold = 0.95,pValueThreshold = 0.05,userLimits = c(-1, 1))


Col22_Qian_Mboi_R <- as.data.table(do.call("rbind",Col22_Qian_Mboi_PCA))
Col_2012_hindiii_R <- as.data.table(do.call("rbind",Col_2012_hindiii_PCA))
Col_Feng_hindiii_R <- as.data.table(do.call("rbind",Col_Feng_hindiii_PCA))
Col_Grob_hindiii_R <- as.data.table(do.call("rbind",Col_Grob_hindiii_PCA))
Col_Chang_dpnii_R <- as.data.table(do.call("rbind",Col_Chang_dpnii_PCA))
Col_Chang_dpnii_2018_R <- as.data.table(do.call("rbind",Col_Chang_dpnii_2018_PCA))




#as.data.table(summary(Matrix(binMatList$Col22_Qian_Mboi, sparse = T))) %>% fwrite(file = "sparse_matrix.tsv", sep = "\t")

full_join(Col22_Qian_Mboi_R, Col_2012_hindiii_R, by = "validBins") %>%
    full_join(Col_Feng_hindiii_R, by = "validBins") %>%
    full_join(Col_Grob_hindiii_R, by = "validBins") %>%
    full_join(Col_Chang_dpnii_R, by = "validBins") %>%
    full_join(Col_Chang_dpnii_2018_R, by = "validBins") %>% as.data.table() -> Col

colnames(Col)[-1] <- c("Col-0 Control",
                       "Col-0 (G. Moissiard)",
                       "Col-0 (S. Feng)",
                       "Col-0 (S. Grob)",
                       "Col-0 (C. Wang)",
                       "Col-0 (W. Zhu)")

Col
## Heatmap of pair-wise PCC of PC1  ----
MM<-as.data.frame(Col)
row.names(MM)<-MM$validBins
MM$validBins<-NULL
MM %>% head()

cor(MM, use = "complete.obs") %>% pheatmap(
    cluster_rows = T,
    cluster_cols = T,
    color = viridis_pal(option = "D")(9),
    border_color = "white",
    fontsize = 15,
    display_numbers = T,
    number_color = "white",
    filename = "PCCofPC1_OfPublishedAndInhouseHiC.pdf",
    height = 10,cellwidth = 60,cellheight = 60,
    width = 10
)

cor(MM, use = "complete.obs") %>% pheatmap(
    cluster_rows = F,
    cluster_cols = F,
    color = viridis_pal(option = "D")(9),
    border_color = "white",
    fontsize = 15,
    display_numbers = T,
    number_color = "white",
    filename = "UnclisteredPCCofPC1_OfPublishedAndInhouseHiC.pdf",
    height = 10,cellwidth = 60,cellheight = 60,
    width = 10
)

## ----

mCol<-melt(Col,id.vars = "validBins",variable.name = "Source",value.name = "Eigenvalue")

mCol$enzyme<-""
mCol[Source=="Col-0 Control"]$enzyme<-"MboI"
mCol[Source=="Col-0 (G. Moissiard)"]$enzyme<-"HindIII"
mCol[Source=="Col-0 (S. Feng)"]$enzyme<-"HindIII"
mCol[Source=="Col-0 (S. Grob)"]$enzyme<-"HindIII"
mCol[Source=="Col-0 (C. Wang)"]$enzyme<-"DpnII"
mCol[Source=="Col-0 (W. Zhu)"]$enzyme<-"DpnII"

mCol$Chr<-""
mCol[validBins<=305]$Chr<-"Chr1 right arm"
mCol[validBins >= 803 & validBins<=923]$Chr<-"Chr4 right arm"
mCol[validBins >= 1084 & validBins<=1193]$Chr<-"Chr5 right arm"


BED<-fread("/data1/linhua/QIANLAB/PROJECT/hic/FPC/BIN100K_Col_Matrix/BU/Col22_total_sample_100000_abs.bed")
colnames(BED)<-c("chr","start","end","validBins")
inner_join(BED,mCol,by = "validBins") %>% as.data.table() -> mCol

pdf(file = "Compare_FPC_All_Col-0.pdf",
    height = 12,
    width = 8)

ggplot(mCol, aes(x = start, y = Eigenvalue, color = Source)) +
    geom_line(aes(linetype = enzyme)) +
    #geom_point(size = 0.5,shape = 21,color = "black",fill = "white") +
    scale_linetype_manual(values = c("longdash", "dotted", "solid")) +
    scale_color_manual(values = c(
        "#EE0000FF",
        "#3B4992FF",
        "#008B45FF",
        "#5F559BFF" ,
        "#808180FF",
        "blue"
    )) +
    scale_y_continuous(limits = c(-0.2, 0.2)) +
    scale_x_log10(name = "Genomic distance (bp)",
                  #breaks = trans_breaks("log10", function(x) 10 ^ x),
                  breaks = base_breaks(n = 4),
                  labels = f2si) -> P

facet(P,
      facet.by = c("Chr"),
      scales = "free_x",
      nrow = 3)
dev.off() 

## Annotation based ----
annotation <- f.read.annotation(file.path("/data1/linhua/QIANLAB/PROJECT/hic/MWSchmid_test/At_tutorial_files", "At_fragments_100kb_annotated.txt"), binSize)
annotationFromFragments <- f.read.annotation.via.fragment.annotation(file.path("/data1/linhua/QIANLAB/PROJECT/hic/MWSchmid_test/At_tutorial_files", "At_fragments_hindIII_annotated.txt"), binSize)

## Correlation analysis with 100K BIN----

Cor_Enrich_Col22_Qian_Mboi<-f.principle.component.analysis.and.features(dataMatrix = binMatList$Col22_Qian_Mboi,outfilePrefix = "Cor_Enrich_Col22_Qian_Mboi",binSize = 100000,rDir = "./",annotation = annotation,regionTable = definedGenomicRegionsArms,simplifiedNames = list(),filterZero = TRUE,filterThreshold = 0.95,pValueThreshold = 0.05,userLimits = c(-1, 1))
Cor_Enrich_Col_2012_hindiii<-f.principle.component.analysis.and.features(dataMatrix = binMatList$Col_2012_hindiii,outfilePrefix = "Cor_Enrich_Col_2012_hindiii",binSize = 100000,rDir = "./",annotation = annotation,regionTable = definedGenomicRegionsArms,simplifiedNames = list(),filterZero = TRUE,filterThreshold = 0.95,pValueThreshold = 0.05,userLimits = c(-1, 1))
Cor_Enrich_Col_Feng_hindiii<-f.principle.component.analysis.and.features(dataMatrix = binMatList$Col_Feng_hindiii,outfilePrefix = "Cor_Enrich_Col_Feng_hindiii",binSize = 100000,rDir = "./",annotation = annotation,regionTable = definedGenomicRegionsArms,simplifiedNames = list(),filterZero = TRUE,filterThreshold = 0.95,pValueThreshold = 0.05,userLimits = c(-1, 1))
Cor_Enrich_Col_Grob_hindiii<-f.principle.component.analysis.and.features(dataMatrix = binMatList$Col_Grob_hindiii,outfilePrefix = "Cor_Enrich_Col_Grob_hindiii",binSize = 100000,rDir = "./",annotation = annotation,regionTable = definedGenomicRegionsArms,simplifiedNames = list(),filterZero = TRUE,filterThreshold = 0.95,pValueThreshold = 0.05,userLimits = c(-1, 1))
Cor_Enrich_Col_Weigel_dpnii<-f.principle.component.analysis.and.features(dataMatrix = binMatList$Col_Chang_dpnii,outfilePrefix = "Cor_Enrich_Col_Weigel_dpnii",binSize = 100000,rDir = "./",annotation = annotation,regionTable = definedGenomicRegionsArms,simplifiedNames = list(),filterZero = TRUE,filterThreshold = 0.95,pValueThreshold = 0.05,userLimits = c(-1, 1))
Cor_Enrich_Col_Liu_dpnii<-f.principle.component.analysis.and.features(dataMatrix = binMatList$Col_Chang_dpnii_2018,outfilePrefix = "Cor_Enrich_Col_Liu_dpnii",binSize = 100000,rDir = "./",annotation = annotation,regionTable = definedGenomicRegionsArms,simplifiedNames = list(),filterZero = TRUE,filterThreshold = 0.95,pValueThreshold = 0.05,userLimits = c(-1, 1))


Cor_Enrich_Col22_Qian_Mboi$estimate[Cor_Enrich_Col22_Qian_Mboi$estimatePvalue>0.05]<-NA
Cor_Enrich_Col_2012_hindiii$estimate[Cor_Enrich_Col_2012_hindiii$estimatePvalue>0.05]<-NA
Cor_Enrich_Col_Feng_hindiii$estimate[Cor_Enrich_Col_Feng_hindiii$estimatePvalue>0.05]<-NA
Cor_Enrich_Col_Grob_hindiii$estimate[Cor_Enrich_Col_Grob_hindiii$estimatePvalue>0.05]<-NA
Cor_Enrich_Col_Weigel_dpnii$estimate[Cor_Enrich_Col_Weigel_dpnii$estimatePvalue>0.05]<-NA
Cor_Enrich_Col_Liu_dpnii$estimate[Cor_Enrich_Col_Liu_dpnii$estimatePvalue>0.05]<-NA

Chr1 <- data.frame(
    Col22_Qian_Mboi_Chr1 = t(Cor_Enrich_Col22_Qian_Mboi$estimate)[, 1],
    Col_2012_hindiii_Chr1 = t(Cor_Enrich_Col_2012_hindiii$estimate)[, 1],
    Col_Feng_hindiii_Chr1 = t(Cor_Enrich_Col_Feng_hindiii$estimate)[, 1],
    Col_Grob_hindiii_Chr1 = t(Cor_Enrich_Col_Grob_hindiii$estimate)[, 1],
    Col_Weigel_dpnii_Chr1 = t(Cor_Enrich_Col_Weigel_dpnii$estimate)[, 1],
    Col_Liu_dpnii_Chr1 = t(Cor_Enrich_Col_Liu_dpnii$estimate)[, 1]
)
Chr4 <- data.frame(
    Col22_Qian_Mboi_Chr4 = t(Cor_Enrich_Col22_Qian_Mboi$estimate)[, 2],
    Col_2012_hindiii_Chr4 = t(Cor_Enrich_Col_2012_hindiii$estimate)[, 2],
    Col_Feng_hindiii_Chr4 = t(Cor_Enrich_Col_Feng_hindiii$estimate)[, 2],
    Col_Grob_hindiii_Chr4 = t(Cor_Enrich_Col_Grob_hindiii$estimate)[, 2],
    Col_Weigel_dpnii_Chr4 = t(Cor_Enrich_Col_Weigel_dpnii$estimate)[, 2],
    Col_Liu_dpnii_Chr4 = t(Cor_Enrich_Col_Liu_dpnii$estimate)[, 2]
)
Chr5 <- data.frame(
    Col22_Qian_Mboi_Chr5 = t(Cor_Enrich_Col22_Qian_Mboi$estimate)[, 3],
    Col_2012_hindiii_Chr5 = t(Cor_Enrich_Col_2012_hindiii$estimate)[, 3],
    Col_Feng_hindiii_Chr5 = t(Cor_Enrich_Col_Feng_hindiii$estimate)[, 3],
    Col_Grob_hindiii_Chr5 = t(Cor_Enrich_Col_Grob_hindiii$estimate)[, 3],
    Col_Weigel_dpnii_Chr5 = t(Cor_Enrich_Col_Weigel_dpnii$estimate)[, 3],
    Col_Liu_dpnii_Chr5 = t(Cor_Enrich_Col_Liu_dpnii$estimate)[, 3]
)
Chr<-cbind(Chr1,Chr4,Chr5)

colnames(Chr)<-c("Col-0 Control Chr1","Col-0 (G. Moissiard) Chr1","Col-0 (S. Feng) Chr1","Col-0 (S. Grob) Chr1","Col-0 (C. Wang) Chr1","Col-0 (W. Zhu) Chr1","Col-0 Control Chr4","Col-0 (G. Moissiard) Chr4","Col-0 (S. Feng) Chr4","Col-0 (S. Grob) Chr4","Col-0 (C. Wang) Chr4","Col-0 (W. Zhu) Chr4","Col-0 Control Chr5","Col-0 (G. Moissiard) Chr5","Col-0 (S. Feng) Chr5","Col-0 (S. Grob) Chr5","Col-0 (C. Wang) Chr5","Col-0 (W. Zhu) Chr5")

c("ann_gene",  "sum_GSM701934_Transcriptome",  "ann_transposable_element",  "ann_si_gene",  "den_GSM701923_H3K4me2",  "den_GSM701924_H3K4me3",  "den_GSM701925_H3K9Ac",  "den_GSM701926_H3K9me2",  "den_GSM701927_H3K18Ac",  "den_GSM701928_H3K27me1",  "den_GSM701929_H3K27me3",  "den_GSM701930_H3K36me2",  "den_GSM701931_H3K36me3",  "den_GSM701932_H3",  "den_CG_rep1",  "den_CG_rep2",  "den_CHG_rep1",  "den_CHG_rep2",  "den_CHH_rep1",  "den_CHH_rep2",  "sum_artificialGenomeReads_sorted")-> SS
Chr[c(rownames(Chr) %in% SS),] -> Chr
rownames(Chr)<-c("gene","transposon","smRNA","transcription","H3K4me2","H3K4me3","H3K9Ac","H3K9me2","H3K18Ac","H3K27me1","H3K27me3","H3K36me2","H3K36me3","H3","CG_rep1","CG_rep2","CHG_rep1","CHG_rep2","CHH_rep1","CHH_rep2","gDNA")


Col_ANN<-data.frame(chromosome=c("Chr1","Chr1","Chr1","Chr1","Chr1","Chr1","Chr4","Chr4","Chr4","Chr4","Chr4","Chr4","Chr5","Chr5","Chr5","Chr5","Chr5","Chr5"))

rownames(Col_ANN)<-colnames(Chr)

pheatmap(
    Chr,
    cluster_cols = F,
    fontsize = 10,
    cellwidth = 12,
    cellheight = 12,border_color = "white",
    cluster_rows = F,
    gaps_col = c(6, 12),
    color = colorRampPalette(rev(brewer.pal(
        n = 9, name =
            "RdBu"
    )))(100),annotation_col = Col_ANN,
    filename = "Correlation_Compare_Diff_Col.pdf",
    height = 6
)

## Enrichment analysis with 100K BIN ----
Cor_Enrich_Col22_Qian_Mboi$enrichment[Cor_Enrich_Col22_Qian_Mboi$enrichmentPvalue>0.05]<-NA
Cor_Enrich_Col_2012_hindiii$enrichment[Cor_Enrich_Col_2012_hindiii$enrichmentPvalue>0.05]<-NA
Cor_Enrich_Col_Feng_hindiii$enrichment[Cor_Enrich_Col_Feng_hindiii$enrichmentPvalue>0.05]<-NA
Cor_Enrich_Col_Grob_hindiii$enrichment[Cor_Enrich_Col_Grob_hindiii$enrichmentPvalue>0.05]<-NA
Cor_Enrich_Col_Weigel_dpnii$enrichment[Cor_Enrich_Col_Weigel_dpnii$enrichmentPvalue>0.05]<-NA
Cor_Enrich_Col_Liu_dpnii$enrichment[Cor_Enrich_Col_Liu_dpnii$enrichmentPvalue>0.05]<-NA

Chr1_E <- data.frame(
    Col22_Qian_Mboi_Chr1 = t(Cor_Enrich_Col22_Qian_Mboi$enrichment)[, 1],
    Col_2012_hindiii_Chr1 = t(Cor_Enrich_Col_2012_hindiii$enrichment)[, 1],
    Col_Feng_hindiii_Chr1 = t(Cor_Enrich_Col_Feng_hindiii$enrichment)[, 1],
    Col_Grob_hindiii_Chr1 = t(Cor_Enrich_Col_Grob_hindiii$enrichment)[, 1],
    Col_Weigel_dpnii_Chr1 = t(Cor_Enrich_Col_Weigel_dpnii$enrichment)[, 1],
    Col_Liu_dpnii_Chr1 = t(Cor_Enrich_Col_Liu_dpnii$enrichment)[, 1]
)
Chr4_E <- data.frame(
    Col22_Qian_Mboi_Chr4 = t(Cor_Enrich_Col22_Qian_Mboi$enrichment)[, 2],
    Col_2012_hindiii_Chr4 = t(Cor_Enrich_Col_2012_hindiii$enrichment)[, 2],
    Col_Feng_hindiii_Chr4 = t(Cor_Enrich_Col_Feng_hindiii$enrichment)[, 2],
    Col_Grob_hindiii_Chr4 = t(Cor_Enrich_Col_Grob_hindiii$enrichment)[, 2],
    Col_Weigel_dpnii_Chr4 = t(Cor_Enrich_Col_Weigel_dpnii$enrichment)[, 2],
    Col_Liu_dpnii_Chr4 = t(Cor_Enrich_Col_Liu_dpnii$enrichment)[, 2]
)
Chr5_E <- data.frame(
    Col22_Qian_Mboi_Chr5 = t(Cor_Enrich_Col22_Qian_Mboi$enrichment)[, 3],
    Col_2012_hindiii_Chr5 = t(Cor_Enrich_Col_2012_hindiii$enrichment)[, 3],
    Col_Feng_hindiii_Chr5 = t(Cor_Enrich_Col_Feng_hindiii$enrichment)[, 3],
    Col_Grob_hindiii_Chr5 = t(Cor_Enrich_Col_Grob_hindiii$enrichment)[, 3],
    Col_Weigel_dpnii_Chr5 = t(Cor_Enrich_Col_Weigel_dpnii$enrichment)[, 3],
    Col_Liu_dpnii_Chr5 = t(Cor_Enrich_Col_Liu_dpnii$enrichment)[, 3]
)

Chr_E<-cbind(Chr1_E,Chr4_E,Chr5_E)
colnames(Chr_E)<-c("Col-0 Control Chr1","Col-0 (G. Moissiard) Chr1","Col-0 (S. Feng) Chr1","Col-0 (S. Grob) Chr1","Col-0 (C. Wang) Chr1","Col-0 (W. Zhu) Chr1","Col-0 Control Chr4","Col-0 (G. Moissiard) Chr4","Col-0 (S. Feng) Chr4","Col-0 (S. Grob) Chr4","Col-0 (C. Wang) Chr4","Col-0 (W. Zhu) Chr4","Col-0 Control Chr5","Col-0 (G. Moissiard) Chr5","Col-0 (S. Feng) Chr5","Col-0 (S. Grob) Chr5","Col-0 (C. Wang) Chr5","Col-0 (W. Zhu) Chr5")

c("ann_gene",  "sum_GSM701934_Transcriptome",  "ann_transposable_element",  "ann_si_gene",  "den_GSM701923_H3K4me2",  "den_GSM701924_H3K4me3",  "den_GSM701925_H3K9Ac",  "den_GSM701926_H3K9me2",  "den_GSM701927_H3K18Ac",  "den_GSM701928_H3K27me1",  "den_GSM701929_H3K27me3",  "den_GSM701930_H3K36me2",  "den_GSM701931_H3K36me3",  "den_GSM701932_H3",  "den_CG_rep1",  "den_CG_rep2",  "den_CHG_rep1",  "den_CHG_rep2",  "den_CHH_rep1",  "den_CHH_rep2",  "sum_artificialGenomeReads_sorted")-> SS
Chr_E[c(rownames(Chr_E) %in% SS),] -> Chr_E
rownames(Chr_E)<-c("gene","transposon","smRNA","transcription","H3K4me2","H3K4me3","H3K9Ac","H3K9me2","H3K18Ac","H3K27me1","H3K27me3","H3K36me2","H3K36me3","H3","CG_rep1","CG_rep2","CHG_rep1","CHG_rep2","CHH_rep1","CHH_rep2","gDNA")


Col_ANN<-data.frame(chromosome=c("Chr1","Chr1","Chr1","Chr1","Chr1","Chr1","Chr4","Chr4","Chr4","Chr4","Chr4","Chr4","Chr5","Chr5","Chr5","Chr5","Chr5","Chr5"))

rownames(Col_ANN)<-colnames(Chr_E)
pheatmap(
    Chr_E,
    cluster_cols = F,
    fontsize = 10,
    cellwidth = 12,
    cellheight = 12,border_color = "white",
    cluster_rows = F,
    gaps_col = c(6, 12),
    color = colorRampPalette(rev(brewer.pal(
        n = 9, name =
            "RdBu"
    )))(100),annotation_col = Col_ANN
    ,
    filename = "Enrichment_Compare_Diff_Col.pdf",
    height = 6
)
## Set up Ploting ----

# Overwrite default draw_colnames in the pheatmap package.
# Thanks to Josh O'Brien at http://stackoverflow.com/questions/15505607
draw_colnames_45 <- function (coln, gaps, ...) {
    coord <- pheatmap:::find_coordinates(length(coln), gaps)
    x     <- coord$coord - 0.5 * coord$size
    res   <- grid::textGrob(
        coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
        vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
    )
    return(res)
}
assignInNamespace(
    x = "draw_colnames",
    value = "draw_colnames_45",
    ns = asNamespace("pheatmap")
)
## Corrrelation analysis for Col (from different sources)----

f.HiC.correlation.matrix(
    dataMatrixList = binMatList,
    rDir = "./",
    outfile = "test_sample",
    corMethod = "pearson",
    summaryFunction = median,
    useOnlyHighVar = T
)

MAT<-read.csv("test_sample_cor_mat.txt",row.names = 1)

colnames(MAT)<-c("Col-0 Control","Col-0 (G. Moissiard)","Col-0 (S. Feng)","Col-0 (S. Grob)","Col-0 (C. Wang)","Col-0 (W. Zhu)")
rownames(MAT)<-c("Col-0 Control","Col-0 (G. Moissiard)","Col-0 (S. Feng)","Col-0 (S. Grob)","Col-0 (C. Wang)","Col-0 (W. Zhu)")

MAT %>% pheatmap(
    cluster_rows = T,
    cluster_cols = T,
    color = viridis_pal(option = "D")(9),
    border_color = "white",
    fontsize = 15,
    display_numbers = T,
    number_color = "white",
    filename = "f.HiC.correlation.matrixOfPublishedAndInhouseHiC.pdf",
    height = 10,cellwidth = 60,cellheight = 60,
    width = 10
)


MAT %>% pheatmap(
    cluster_rows = F,
    cluster_cols = F,
    color = viridis_pal(option = "D")(9),
    border_color = "white",
    fontsize = 15,
    display_numbers = T,
    number_color = "white",
    filename = "UnClusteredf.HiC.correlation.matrixOfPublishedAndInhouseHiC.pdf",
    height = 10,cellwidth = 60,cellheight = 60,
    width = 10
)
