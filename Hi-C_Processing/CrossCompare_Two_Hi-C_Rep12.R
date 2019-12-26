## Import R packages and Define R functions ----
library(pacman)
suppressWarnings(suppressPackageStartupMessages(
    p_load(
        HiCdatR,
        data.table,
        ggplot2,
        viridis,
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
    #strip.background = element_blank(),
    #panel.border = element_blank(),
    #aspect.ratio = 2 / (1 + sqrt(5)),
    aspect.ratio = 1,
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

setwd("/data1/linhua/QIANLAB/PROJECT/hic/FPC/BIN100K_Col_Matrix/Lib12_compare")
pathToTutorial <- "/data1/linhua/QIANLAB/PROJECT/hic/FPC/BIN100K_Col_Matrix/Lib12_compare"
pathToScripts <- "/data1/linhua/QIANLAB/PROJECT/hic/MWSchmid_test/Rscripts/"
f.source.organism.specific.code(file.path(pathToScripts, "HiCdat-A-thaliana-TAIR10.R"))
definedGenomicRegionsArms <- read.table(file.path("/data1/linhua/QIANLAB/PROJECT/hic/MWSchmid_test/At_tutorial_files", "genomicRegionsArms_forPCA.txt"), sep = '\t', header = TRUE, quote = "", stringsAsFactors = FALSE)
binSize <- 20000


sampleList <- list(
    Control_Lib1 = c("Col22_Lib1_sample_20000.matrix"),
    Control_Lib2 = c("Col22_Lib2_sample_20000.matrix"),
    Heat_Lib1 = c("Col37_Lib1_sample_20000.matrix"),
    Heat_Lib2 = c("Col37_Lib2_sample_20000.matrix")
)
binMatList <- f.load.samples(pathToTutorial, sampleList, binSize, 50) ## 5959*5959 MATRIX
#saveRDS(binMatList,file = "Col_binMatList_Lib12.RDS")

Control_Lib1 <-
    f.principle.component.analysis.and.features(
        dataMatrix = binMatList$Control_Lib1,
        outfilePrefix = "Control_Lib1",
        binSize = 20000,
        rDir = "./",
        regionTable = definedGenomicRegionsArms,
        simplifiedNames = list(),
        filterZero = TRUE,
        filterThreshold = 0.95,
        pValueThreshold = 0.05,
        userLimits = c(-1, 1)
    )

Control_Lib2 <-
    f.principle.component.analysis.and.features(
        dataMatrix = binMatList$Control_Lib2,
        outfilePrefix = "Control_Lib2",
        binSize = 20000,
        rDir = "./",
        regionTable = definedGenomicRegionsArms,
        simplifiedNames = list(),
        filterZero = TRUE,
        filterThreshold = 0.95,
        pValueThreshold = 0.05,
        userLimits = c(-1, 1)
    )
Heat_Lib1 <-
    f.principle.component.analysis.and.features(
        dataMatrix = binMatList$Heat_Lib1,
        outfilePrefix = "Heat_Lib1",
        binSize = 20000,
        rDir = "./",
        regionTable = definedGenomicRegionsArms,
        simplifiedNames = list(),
        filterZero = TRUE,
        filterThreshold = 0.95,
        pValueThreshold = 0.05,
        userLimits = c(-1, 1)
    )
Heat_Lib2 <-
    f.principle.component.analysis.and.features(
        dataMatrix = binMatList$Heat_Lib2,
        outfilePrefix = "Heat_Lib2",
        binSize = 20000,
        rDir = "./",
        regionTable = definedGenomicRegionsArms,
        simplifiedNames = list(),
        filterZero = TRUE,
        filterThreshold = 0.95,
        pValueThreshold = 0.05,
        userLimits = c(-1, 1)
    )

Control_Lib1_DT <- as.data.table(do.call("rbind",Control_Lib1))
Control_Lib2_DT <- as.data.table(do.call("rbind",Control_Lib2))
Heat_Lib1_DT <- as.data.table(do.call("rbind",Heat_Lib1))
Heat_Lib2_DT <- as.data.table(do.call("rbind",Heat_Lib2))

full_join(Control_Lib1_DT, Control_Lib2_DT, by = "validBins") %>%
    full_join(Heat_Lib1_DT, by = "validBins") %>%
    full_join(Heat_Lib2_DT, by = "validBins") %>%
    as.data.table() -> Col

colnames(Col)[-1] <- c("Control_Rep1","Control_Rep2","Heat_Rep1","Heat_Rep2")
mCol<-melt(Col,id.vars = "validBins",variable.name = "Source",value.name = "Eigenvalue")
mCol$Chr<-""

BED<-fread("/data1/linhua/QIANLAB/PROJECT/hic/Matrix_Heat/Col22_total_sample/raw/20000/Col22_total_sample_20000_abs.bed")
colnames(BED)<-c("chr","start","end","validBins")
inner_join(BED,mCol,by = "validBins") %>% as.data.table() -> mCol
mCol[chr=="Chr1"]$Chr<-"Chr1 right arm"
mCol[chr=="Chr4"]$Chr<-"Chr4 right arm"
mCol[chr=="Chr5"]$Chr<-"Chr5 right arm"


pdf(file = "Compare_FPC_All_LIB12.pdf",
    height = 12,
    width = 8)
library(RColorBrewer)

ggplot(mCol, aes(x = start, y = Eigenvalue, color = Source)) +
    geom_line(aes(linetype = Source)) +
    #geom_point(size = 0.5,shape = 21,color = "black",fill = "white") +
    scale_linetype_manual(values = c("dotted", "solid","dotted", "solid")) +
    scale_colour_manual(values = brewer.pal(6, "Paired")[c(1,2,5,6)])+
    scale_y_continuous(limits = c(-0.1, 0.1)) +
    scale_x_log10(name = "Genomic distance (bp)",
                  #breaks = trans_breaks("log10", function(x) 10 ^ x),
                  breaks = base_breaks(n = 4),
                  labels = f2si) -> P
facet(P,
      facet.by = c("Chr"),
      scales = "free_x",
      nrow = 3)

dev.off()

## Pair-wise scatterplot comparision of PC1 ----
Get2D_KDEP <- function(Data,T1, T2) {
    ggplot(Data, aes_string(T1, T2)) + 
        geom_point(size = 0.8,color = "black",fill = "white",shape = 21)+
        stat_cor(method = "pearson",color = "black",size = 6,label.x=-0.05,label.y = 0.075) +
        #stat_density_2d(aes(fill = ..level..), geom = "polygon",n = 50) +
        #scale_fill_viridis(option = "A") +
        scale_x_continuous(name = T1,expand = c(0, 0),limits = c(-0.1, 0.1)) +
        scale_y_continuous(name = T2,expand = c(0, 0),limits = c(-0.1, 0.1)) +
        geom_abline(intercept = 0,slope = 1,color = 'blue',size = 0.7)
}

Get2D_KDEP(Data = Col,T1 = "Control_Rep1",T2 = "Control_Rep2") 	%>% ggsave(filename = "Control_Rep1_Control_Rep2.pdf",device = "pdf",height=6,width=6)
Get2D_KDEP(Data = Col,T1 = "Heat_Rep1",T2 = "Heat_Rep2")	%>% ggsave(filename = "Heat_Rep1_Heat_Rep2.pdf",device = "pdf",height=6,width=6)
Get2D_KDEP(Data = Col,T1 = "Control_Rep1",T2 = "Heat_Rep1")	%>% ggsave(filename = "Control_Rep1_Heat_Rep1.pdf",device = "pdf",height=6,width=6)
Get2D_KDEP(Data = Col,T1 = "Control_Rep1",T2 = "Heat_Rep2")	%>% ggsave(filename = "Control_Rep1_Heat_Rep2.pdf",device = "pdf",height=6,width=6)
Get2D_KDEP(Data = Col,T1 = "Control_Rep2",T2 = "Heat_Rep1")	%>% ggsave(filename = "Control_Rep2_Heat_Rep1.pdf",device = "pdf",height=6,width=6)
Get2D_KDEP(Data = Col,T1 = "Control_Rep2",T2 = "Heat_Rep2")	%>% ggsave(filename = "Control_Rep2_Heat_Rep2.pdf",device = "pdf",height=6,width=6)


## Heatmap of pair-wise PCC of PC1  ----
MM<-as.data.frame(Col)
row.names(MM)<-MM$validBins
MM$validBins<-NULL
MM %>% head()

cor(MM, use = "complete.obs") %>% pheatmap(
    cluster_rows = F,
    cluster_cols = F,
    color = viridis_pal(option = "D")(9),
    border_color = "white",
    fontsize = 15,
    display_numbers = T,
    number_color = "white",
    filename = "ControlRep12AndHeatRep12.pdf",
    height = 6,cellwidth = 60,cellheight = 60,
    width = 6
)

# Use f.HiC.correlation.matrix for Hi-C QC ----

f.HiC.correlation.matrix(
    dataMatrixList = binMatList,
    rDir = "./",
    outfile = "Control_Heat_lib12",
    corMethod = "spearman",
    summaryFunction = median,
    useOnlyHighVar = T
)
MAT<-read.csv("Control_Heat_lib12_cor_mat.txt",row.names = 1)

row.names(MAT)<-str_replace_all(row.names(MAT),"Lib","Rep")
base::colnames(MAT)<-str_replace_all(base::colnames(MAT),"Lib","Rep")

library(pheatmap)
pheatmap(
    MAT,
    cluster_rows = F,
    cluster_cols = F,
    color = viridis_pal(option = "D")(9),
    border_color = "white",
    fontsize = 15,
    display_numbers = T,
    number_color = "white",
    filename = "ControlRep12AndHeatRep12.HiC.correlation.pdf",
    height = 6,cellwidth = 60,cellheight = 60,
    width = 6 
)

## Direct heatmap compare for KNOB region ----

library(Matrix)
Control <- binMatList$Control_Lib1
Treat <- binMatList$Control_Lib2
diag(Control) <- NA
diag(Treat) <- NA
Control[lower.tri(Control)] <- 0
Treat[upper.tri(Treat)] <- 0
Control + Treat -> MAT

as.data.table(summary(Matrix(MAT, sparse = T))) -> MAT
source("/data1/linhua/QIANLAB/PROJECT/hic/gg_Tang_Jiaxiang.R")
pdf(file = "Triangle_Heatmap_Compare_Control_lib12.pdf",
    height = 8,
    width = 8)
ggtrapezoid(
    mat = MAT,
    bed = BED,
    type = "full",
    region = "Chr4:1-6000000"
)
dev.off()