## Library R packages ----
library(pacman)
suppressWarnings(suppressPackageStartupMessages(
    p_load(
        ComplexHeatmap,
        DESeq2,
        viridis,
        knitr,
        crayon,
        data.table,
        dendsort,
        dplyr,
        edgeR,
        extrafont,
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
        ggrepel
    )
))
suppressMessages(loadfonts())

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

## Aim: audo DEG analysis for multiple sample (with >3 replicates)

## Import quant value for each sample ----

Args<-commandArgs(trailingOnly = T)

if (length(Args)==0) {
    stop("At least 1 argument must be supplied: Rscriipt script.R Ctrl1 Ctrl2 Mut1 Mut2", call.=FALSE)
}

#setwd("/data1/linhua/QIANLAB/PROJECT/ddccm/RNA-seq")


Args[1]<-"wt_RNA_rep1_STARAligned.sortedByCoord.out.txt"
Args[2]<-"wt_RNA_rep2_STARAligned.sortedByCoord.out.txt"
Args[3]<-"wt_RNA_rep3_STARAligned.sortedByCoord.out.txt"
Args[4]<-"wt_RNA_rep4_STARAligned.sortedByCoord.out.txt"
Args[5]<-"h1_RNA_rep1_STARAligned.sortedByCoord.out.txt"
Args[6]<-"h1_RNA_rep2_STARAligned.sortedByCoord.out.txt"
Args[7]<-"h1_RNA_rep3_STARAligned.sortedByCoord.out.txt"
Args[8]<-"h1_RNA_rep4_STARAligned.sortedByCoord.out.txt"
Args[9]<-"h1_Daniel"

Ctrl1<-Args[1]
Ctrl2<-Args[2]
Ctrl3<-Args[3]
Ctrl4<-Args[4]
Mut1<-Args[5]
Mut2<-Args[6]
Mut3<-Args[7]
Mut4<-Args[8]
OutKey<-Args[9]

ImportFC <- function(filepath) {
    heat_fc<-fread(filepath,select = c(1,7),skip = 1)
    str_replace_all(colnames(heat_fc),
                    pattern = "uniqmap_|_tophat_align_sortbyname.bam",
                    replacement = "") %>% str_replace_all(pattern = "Geneid", "ID") -> colnames(heat_fc)
    return(heat_fc)
}

Ctrl1_RawCounts<-ImportFC(Ctrl1)
Ctrl2_RawCounts<-ImportFC(Ctrl2)
Ctrl3_RawCounts<-ImportFC(Ctrl3)
Ctrl4_RawCounts<-ImportFC(Ctrl4)
Mut1_RawCounts<-ImportFC(Mut1)
Mut2_RawCounts<-ImportFC(Mut2)
Mut3_RawCounts<-ImportFC(Mut3)
Mut4_RawCounts<-ImportFC(Mut4)

full_join(Ctrl1_RawCounts,Ctrl2_RawCounts,by="ID") %>% 
    full_join(Ctrl3_RawCounts,by="ID") %>% 
    full_join(Ctrl4_RawCounts,by="ID") %>% 
    full_join(Mut1_RawCounts,by="ID") %>% 
    full_join(Mut2_RawCounts,by="ID") %>% 
    full_join(Mut3_RawCounts,by="ID") %>% 
    full_join(Mut4_RawCounts,by="ID") %>% 
    as.data.table() -> heat_fc

back<-heat_fc

colnames(heat_fc)<-c("ID","Ctrl1","Ctrl2","Ctrl3","Ctrl4","Mut1","Mut2","Mut3","Mut4")

## Convert raw rounts into data.frame format (for DEseq2 )----
heat_fc_df<-as.data.frame(heat_fc)
rownames(heat_fc_df)<-heat_fc_df$ID
heat_fc_df$ID<-NULL
#fwrite(heat_fc_df,file = "heat_fc_raw_counts.csv")

## featurecounts reads distribution analysis ----
Ath_GENE_TE_Clean <- fread(file = "/data1/linhua/QIANLAB/PROJECT/te/Ath_GENE_TE_Clean_Pseudo.csv")

full_join(heat_fc,Ath_GENE_TE_Clean,by="ID") %>%
    as.data.table() -> ANN

melt(ANN,id.vars = c("ID","group","locus_type")) -> MANN

MANN %>%
    group_by(group,variable) %>%
    summarise(NO_mapped_reads=sum(value)) %>%
    group_by(variable) %>%
    mutate(fraction=NO_mapped_reads/sum(NO_mapped_reads)*100) -> PP

#PP %>% kable(digits = 3,format = "rst")

PP %>% ggbarplot(x = "variable",
                 y = "fraction",
                 group = "group",width = 0.7,
                 fill = "group") %>%
    ggpar(
        ylim = c(0, round(max(PP$fraction[PP$fraction<10]))+0.5),
        x.text.angle = 45,
        xlab = "",
        ylab = "Fraction of reads distribution (%)"
    ) + border()+
    theme(aspect.ratio = 2) -> Reads_Distribution_plot


## Import into DEseq2 env ----

coldata<-data.frame(treat=c("Control","Control","Control","Control","Mutant","Mutant","Mutant","Mutant"))

rownames(coldata)<-colnames(heat_fc_df)

F1<-heat_fc_df
F1[rowSums(F1)>=40,] -> F1 ## Filter

#log2(rowSums(heat_fc_df)+1) %>% gghistogram(bins = 30) + geom_vline(xintercept = c(log2(20),log2(1),log2(50) ))

if (all(rownames(coldata) %in% colnames(F1))) {
    Unit_dds <- DESeqDataSetFromMatrix(
        countData = round(F1),
        colData = coldata,
        design = ~ treat
    )
    #Unit_dds <- Unit_dds[ rowSums(counts(Unit_dds)) > ncol(F1)-1, ] ## Over 1 each sample
}

rld <- vst(Unit_dds, blind=FALSE)
rlog_EX<-as.data.frame(assay(rld))
#fwrite(rlog_EX,file = "ssRNA_Seq_9_samples_VST_reads_counts.csv",row.names = T)

## PCA analysis of top varied 1000 genes ----
rlog_EX<-as.data.frame(assay(rld))

rlog_EX<-rlog_EX[order(genefilter::rowVars(rlog_EX),decreasing = T),]

# transpose the data
data <- t(as.matrix(rlog_EX[1:1000,])) ##
data.pca <- fast.prcomp(data,scale=T,center = T)
a <- summary(data.pca)
tmp <- a[4]$importance
pro1 <- as.numeric(sprintf("%.3f",tmp[2,1]))*100
pro2 <- as.numeric(sprintf("%.3f",tmp[2,2]))*100
pro3 <- as.numeric(sprintf("%.3f",tmp[2,3]))*100

pc = as.data.frame(a$x)
rownames(pc)-> pc$names
pc$group <- coldata$treat
pc$color<-"black"
pc[pc$group== "Control",]$color<-pal_aaas(alpha = 1)(9)[1]
pc[pc$group== "Mutant",]$color<-pal_aaas(alpha = 1)(9)[2]


## 3D PCA PLOT ----

pdf(file =paste(OutKey,"ssRNA_Seq_PCA_3D_genes_DEseq2.pdf",sep = ""),height = 7,width = 8)
with(pc, {
    s3d <-
        scatterplot3d(
            pc$PC1,
            pc$PC2,
            pc$PC3,
            angle = 55,
            color = pc$color,
            pch = 19,
            scale.y = 1.25,
            cex.symbols = 2.0,
            xlab = paste("PC1 (", pro1, "%)", sep = ""),
            ylab = paste("PC2 (", pro2, "%)", sep = ""),
            zlab = paste("PC3 (", pro3, "%)", sep = "")
        )
    s3d.coords <- s3d$xyz.convert(pc$PC1, pc$PC2, pc$PC3)
    text(s3d.coords$x,s3d.coords$y,labels = str_replace_all(str_replace_all(row.names(pc), "IAA_", ""), "_h", "h"),pos = 4,cex = 1)
    legend("bottomright",inset = 0,bty = "n",cex = 1.2,title = "Group",str_c(unique(pc$group), ""),fill = pal_aaas(alpha = 1)(5))
})
dev.off()

## 2D PCA plot (ggplot2) ----
pdf(file = paste(OutKey,"polyA_RNA_Seq_PCA_2D_genes_DEseq2.pdf",sep = ""),height = 7,width = 8)
ggplot(pc, aes(PC1, PC2)) +
    geom_point(
        aes(fill = group),
        shape = 21,
        colour = "black",
        size = 6,
        stroke = 0.8,
        alpha = 1
    ) +
    scale_fill_ucscgb() +
    theme_pubr() ->P
ggpar(
    P,
    xlab = paste("PC1 (", pro1, "%)", sep = ""),
    ylab = paste("PC2 (", pro2, "%)", sep = ""),
    legend = "top",
    #font.family = "Arial",
    font.x = c("bold", 18),
    font.y = c("bold", 18),
    font.tickslab = c("bold", 15, "black"),
    font.legend = c("bold", 15)
) + border() + theme(aspect.ratio = 2 / (1 + sqrt(5)))# golden ratio portrait
dev.off()

## DESeq2 (Pair-Wise analysis) in batch mode ----

Unit_dds<-DESeq(Unit_dds)

H_VS_C<-as.data.table(as.data.frame(results(Unit_dds,c("treat","Mutant","Control"))),keep.rownames = "ID")[,"direction":="NoChange"][]

temp<-H_VS_C
temp[padj < 0.05 & log2FoldChange>1]$direction<-"UP"
temp[padj < 0.05 & log2FoldChange<(-1)]$direction<-"DOWN"

temp[str_detect(ID,"TE")] %>% group_by(direction) %>% summarise(count=n()) %>% as.data.table() %$% .[,type:="TE"] -> polyA_DE_TE_No
temp[!str_detect(ID,"TE")] %>% group_by(direction) %>% summarise(count=n())%>% as.data.table() %$% .[,type:="PCG"] -> polyA_DE_GENE_No

MM<-rbind(polyA_DE_GENE_No,polyA_DE_TE_No)

ggbarplot(
    data = MM[direction != "NoChange"],
    x = "type",
    y = "count",
    color = "black",
    fill = "direction",
    palette = "aaas",
    label = "count",
    width = 0.5,
    position = position_fill(),
    size = 0.5,
    lab.size = 6,
    lab.col = "white",
    lab.pos = "in"
) %>% ggpar(
    legend = "right",
    xlab = "",
    ylab = "Fraction and number\n of TE and PCG",
    font.x = c(15),
    font.y = c(15),
    font.tickslab = c(12),
    font.legend = c(12),
    legend.title = ""
) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(aspect.ratio = 2, axis.line.x = element_blank()) -> TE_PCG_Frac_plot



## Prepare CPM table from raw counts ----
ExCPM<-as.data.table(cpm(F1),keep.rownames = "ID")

ExCPM$Control.mean.cpm<-rowMeans(ExCPM[,2:5])
ExCPM$Mutant.mean.cpm<-rowMeans(ExCPM[,6:9])

left_join(ExCPM,
          heat_fc,
          by = "ID",
          suffix = c(".cpm", ".rawcount")) %>%
    as.data.table() -> AA

## Table making from CPM and DEG analysis into a big table ----
AA %>% full_join(temp,by="ID") %>% as.data.table() %$% .[,`:=`(Control.mean.cpm.log2=log2(Control.mean.cpm+1),                                                             Mutant.mean.cpm.log2=log2(Mutant.mean.cpm+1)) ,by="ID"][] -> IIIII

#fwrite(IIIII,file = "heat_polyA-RNA-seq-table.csv")

## Annotation of PCG  ----
#GFF<-as.data.table(import("/data1/linhua/QIANLAB/Araport_TOPHAT/Araport11_GFF3_genes_transposons.201606.gff"))
#GFF[type=='gene' | type=="transposable_element_gene" |type=="pseudogene"][,c("type","locus_type","ID","Note"#,"symbol","Alias","full_name","seqnames","start","end","width","strand")] -> LLL
#
#LLL$Alias %>% sapply(paste, collapse = "|",simplify = T) -> LLL$Alias
#LLL$Note %>% sapply(paste, collapse = "|",simplify = T) -> LLL$Note
#fwrite(LLL,"Ath_PCG_Annotation.csv")
LLL<-fread("/data1/linhua/QIANLAB/PROJECT/ddccm/RNA-seq/Ath_PCG_Annotation.csv")
## Calculating RPKM for genes ----
LENGTH <- fread(file = "/data1/linhua/QIANLAB/PROJECT/te/Ath_New_GENE_TE_Length.csv")

left_join(IIIII,LENGTH,by="ID") %>% as.data.table() -> IIIII_ANN

rpkm(IIIII_ANN[,c(colnames(IIIII_ANN)[str_detect(colnames(IIIII_ANN),"count")]),with=F],IIIII_ANN$length) %>% as.data.table() -> RPKM
colnames(RPKM)<-str_replace_all(colnames(RPKM),"rawcount","rpkm")

cbind(IIIII_ANN,RPKM) -> IIIII_RPKM

IIIII_RPKM[, "Ctrl.mean.rpkm":= mean(`Ctrl1.rpkm`+`Ctrl2.rpkm`+`Ctrl3.rpkm`+`Ctrl4.rpkm`),by="ID"][, "Mut.mean.rpkm":= mean(`Mut1.rpkm`+`Mut2.rpkm`+`Mut3.rpkm`+`Mut4.rpkm`),by="ID"][, "Ctrl.mean.rpkm.log2":= log2(Ctrl.mean.rpkm+1),by="ID"][, "Mut.mean.rpkm.log2":= log2(Mut.mean.rpkm+1),by="ID"][]

#IIIII_ANN %>% fwrite(file = "Heat_Stress_polyA-RNA-Seq_All_Ann_Results.csv")


## Combind Genes Expression level, DEG and Gene annotation and scatterplot for Genes ----
inner_join(IIIII_RPKM[!str_detect(ID,"TE")],LLL,by = "ID") %>% as.data.table() -> GENE_ANN


ORDER <- c("ID","direction","Control.mean.cpm","Mutant.mean.cpm","Ctrl.mean.rpkm","Mut.mean.rpkm","log2FoldChange","group","length","type","locus_type","Note","symbol","Alias","full_name","baseMean","lfcSE","stat","pvalue","padj","Control.mean.cpm.log2","Mutant.mean.cpm.log2","Ctrl.mean.rpkm.log2","Mut.mean.rpkm.log2",
           "Ctrl1.rawcount",
           "Ctrl2.rawcount",
           "Ctrl3.rawcount",
           "Ctrl4.rawcount",
           "Mut1.rawcount",
           "Mut2.rawcount",
           "Mut3.rawcount",
           "Mut4.rawcount",
           "Ctrl1.cpm",
           "Ctrl2.cpm",
           "Ctrl3.cpm",
           "Ctrl4.cpm",
           "Mut1.cpm",
           "Mut2.cpm",
           "Mut3.cpm",
           "Mut4.cpm",
           "Ctrl1.rpkm",
           "Ctrl2.rpkm",
           "Ctrl3.rpkm",
           "Ctrl4.rpkm",
           "Mut1.rpkm",
           "Mut2.rpkm",
           "Mut3.rpkm",
           "Mut4.rpkm",
           "seqnames","start","end","width","strand")

GENE_ANN[order(ID,decreasing = F)][, ORDER , with = F] -> GENE_ANN

P <-ggplot(GENE_ANN, aes(x = Ctrl.mean.rpkm.log2 , y = Mut.mean.rpkm.log2)) +
    scale_x_continuous(name = "Control (log2(rpkm))") +
    scale_y_continuous(name = "Mutant (log2(rpkm))")
#limits = c(0, 15) 
#limits = c(0, 15) 
L1 <- P + stat_density_2d(aes(fill = ..level..), geom = "polygon") +
    scale_fill_viridis(option = "D")

L2 <- geom_point(data = GENE_ANN[direction=="UP"],aes(x = Ctrl.mean.rpkm.log2, y = Mut.mean.rpkm.log2),color = "#BE2442",alpha = 0.4,size = 2)

L3 <- geom_point(data = GENE_ANN[direction=="DOWN"],aes(x = Ctrl.mean.rpkm.log2, y = Mut.mean.rpkm.log2),color = "#324743",alpha = 0.4,size = 2)

L1 +L2 + L3 +
    geom_abline(slope = 1, intercept = 0,lty="solid",color = "blue",size = 0.7) +
    border() +
    theme(aspect.ratio = 1) -> GENE_Scatterplot




## TE ANN ----

TE_Annotation <- fread("/data1/linhua/QIANLAB/PROJECT/te/TE_Annotation-2018-11-17.csv")
TE_Annotation[,length:=NULL]
inner_join(IIIII_RPKM[str_detect(ID, "TE")] , TE_Annotation, by="ID") %>% as.data.table() -> TE_RNA


TE_Col_Order <- c("ID","direction","Control.mean.cpm","Mutant.mean.cpm","Ctrl.mean.rpkm","Mut.mean.rpkm","log2FoldChange","group","length","Position","Ol_TEGs","Chromatin","Transposon_Family","Transposon_Super_Family","Class","score.GC","score.H3K9me2","baseMean","lfcSE","stat","pvalue","padj",
                  "Ctrl1.rawcount",
                  "Ctrl2.rawcount",
                  "Ctrl3.rawcount",
                  "Ctrl4.rawcount",
                  "Mut1.rawcount",
                  "Mut2.rawcount",
                  "Mut3.rawcount",
                  "Mut4.rawcount",
                  "Ctrl1.cpm",
                  "Ctrl2.cpm",
                  "Ctrl3.cpm",
                  "Ctrl4.cpm",
                  "Mut1.cpm",
                  "Mut2.cpm",
                  "Mut3.cpm",
                  "Mut4.cpm",
                  "Ctrl1.rpkm",
                  "Ctrl2.rpkm",
                  "Ctrl3.rpkm",
                  "Ctrl4.rpkm",
                  "Mut1.rpkm",
                  "Mut2.rpkm",
                  "Mut3.rpkm",
                  "Mut4.rpkm",
                  "Control.mean.cpm.log2","Mutant.mean.cpm.log2","Ctrl.mean.rpkm.log2","Mut.mean.rpkm.log2","seqnames","start","end")

TE_RNA[,TE_Col_Order,with=F] -> TE_RNA

#fwrite(TE_RNA,file = "Heat_Stress_polyA-RNA-Seq_DEG_TEs_Ann_Results.csv")

#fwrite(TE_RNA[Heat_vs_Control=="UP"],file = "heat_polyA_UP_TE_Ann.csv")
#fwrite(TE_RNA[Heat_vs_Control=="DOWN"],file = "heat_polyA_DOWN_TE_Ann.csv")



## TE up and down regulated scatterplot (Refined version) ----
ggplot(TE_RNA, aes(x = Control.mean.cpm.log2, y = Mutant.mean.cpm.log2,label="ONSEN")) +
    geom_point(
        data = TE_RNA[direction == "NoChange"],
        aes(x = Control.mean.cpm.log2, y = Mutant.mean.cpm.log2),
        color = "black",
        alpha = 0.7,
        size = 3,
        fill = 'lightgrey',
        shape = 21
    ) +
    geom_point(
        data = TE_RNA[direction == "UP"],
        aes(x = Control.mean.cpm.log2, y = Mutant.mean.cpm.log2),
        color = "black",
        alpha = 0.6,
        size = 3,
        fill = 'red',
        shape = 21
    ) +
    geom_point(
        data = TE_RNA[direction == "DOWN"],
        aes(x = Control.mean.cpm.log2, y = Mutant.mean.cpm.log2),
        color = "black",
        alpha = 0.6,
        size = 3,
        fill = 'blue',
        shape = 21
    ) +
    scale_x_continuous(limits = c(0, 10), name = "Control (log2(cpm))") +
    scale_y_continuous(limits = c(0, 10), name = "Mutant (log2(cpm))") +
    geom_abline(intercept = 0, color = "purple",size=1) +
    theme(aspect.ratio = 1) + border() -> ScatterPlot_TE_Pair
#geom_abline(slope = 1,intercept = 1, color = "purple",size=1)+
#geom_abline(slope = 1,intercept = -1, color = "purple",size=1)

#ggsave(filename = "TE_Heat_VS_Control_polyA_ONSEN_Green_label.pdf",height = 6,width = 6)

## rbind Total and Up TE for compare ----

TE_Annotation <- fread("/data1/linhua/QIANLAB/PROJECT/te/TE_Annotation-2018-11-17.csv")
TE_RNA[direction=="UP"][,intersect(colnames(TE_Annotation),colnames(TE_RNA)),with=F][,"group":="Up_TE"]-> TE_ANN
TE_Annotation[,intersect(colnames(TE_Annotation),colnames(TE_RNA)),with=F][,"group":="Total"] -> Total_ANN
rbind(TE_ANN,Total_ANN) -> MM
MM$group<-factor(MM$group)


## Set UP ploting parameters ----
library(pacman)
suppressWarnings(suppressPackageStartupMessages(
    p_load(
        ggpubr,
        ggplot2,
        ggsci,
        stringr,
        data.table,
        hashmap,
        crayon,
        extrafont
    )
))
suppressMessages(loadfonts())

source("/data1/linhua/QIANLAB/R_Functions/theme_linhua.R")

theme_set(theme_linhua(base_size = 18,base_family = "",legend = "none",border = F))
theme_update(panel.grid.minor.x = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.grid.major.x = element_blank(),
             panel.grid.major.y = element_blank(),
             axis.line.x = element_blank(),
             axis.title.x = element_text(size = 15),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             panel.spacing.x = unit(0.5, "lines"),
             panel.spacing.y = unit(0.5, "lines"),
             strip.background = element_blank(),
             panel.border = element_blank(),
             axis.line.y = element_line(colour = "black"),
             aspect.ratio = 2
)

## UP TE length, GC contents and H3K9me2 enrichment analysis ----

MM$score.H3K9me2.log10<-log10(MM$score.H3K9me2+1)
ggplot(MM[!is.na(score.H3K9me2.log10)], aes(x = group, y = score.H3K9me2.log10 , fill = group)) +
    stat_boxplot(geom = 'errorbar', width = 0.4) +
    geom_boxplot(outlier.colour = "black",width = 0.6,outlier.size = 1,notch = T,notchwidth = 0,outlier.shape = 1) +
    stat_summary(fun.y = median,geom = "line",aes(group = 1),size = 1,color = "deepskyblue2") +
    stat_summary(fun.y = median,geom = "point",shape = 21,size = 3,fill = "yellow",alpha = 0.8) +
    scale_fill_aaas(alpha = 1) +
    labs(x = "", y = "H3K9me2 score (log10)") +
    scale_y_continuous(expand = c(0, 0)) +
    stat_compare_means(label.x =1.5,label.y = 1.5,aes(label = paste0("p ", ..p.format..)),size=5)+border() -> H3K9me2_plot

ggplot(MM[!is.na(score.GC)], aes(x = group, y = score.GC , fill = group)) +
    stat_boxplot(geom = 'errorbar', width = 0.4) +
    geom_boxplot(outlier.colour = "black",width = 0.6,outlier.size = 1,notch = T,notchwidth = 0,outlier.shape = 1) +
    stat_summary(fun.y = median,geom = "line",aes(group = 1),size = 1,color = "deepskyblue2") +
    stat_summary(fun.y = median,geom = "point",shape = 21,size = 3,fill = "yellow",alpha = 0.8) +
    scale_fill_aaas(alpha = 1) +
    labs(x = "", y = "GC contents") +
    scale_y_continuous(expand = c(0, 0)) +
    stat_compare_means(label.x =1.5,label.y = 0.6,aes(label = paste0("p ", ..p.format..)),size=5)+border() -> GC_plot

MM$length_log10<-log10(MM$length)
ggplot(MM[!is.na(length_log10)], aes(x = group, y = length_log10 , fill = group)) +
    stat_boxplot(geom = 'errorbar', width = 0.4) +
    geom_boxplot(outlier.colour = "black",width = 0.6,outlier.size = 1,notch = T,notchwidth = 0,outlier.shape = 1) +
    stat_summary(fun.y = median,geom = "line",aes(group = 1),size = 1,color = "deepskyblue2") +
    stat_summary(fun.y = median,geom = "point",shape = 21,size = 3,fill = "yellow",alpha = 0.8) +
    scale_fill_aaas(alpha = 1) +
    labs(x = "", y = "length (log10)") +
    scale_y_continuous(expand = c(0, 0)) +
    stat_compare_means(label.x =1.5,label.y = 4,aes(label = paste0("p ", ..p.format..)),size=5) + border() -> Length_plot

## Output----

H3K9me2_plot %>%  ggsave(filename = paste(OutKey,"_H3K9me2_log10_TE_Enrichment_Boxplot.pdf",sep=""),height = 8,width = 4)
GC_plot %>%  ggsave(filename = paste(OutKey,"_GC_contents_TE_Enrichment_Boxplot.pdf",sep=""),height = 8,width = 4)
Length_plot %>%  ggsave(filename = paste(OutKey,"_Length_log10_TE_Enrichment_Boxplot.pdf",sep=""),height = 8,width = 4)



TE_PCG_Frac_plot %>%  ggsave(filename = paste(OutKey,"_TE_PCG_Frac_plot.pdf",sep=""),height = 8,width = 4)
Reads_Distribution_plot %>%  ggsave(filename = paste(OutKey,"_reads_distribution.pdf",sep=""),height = 8,width = 4)



PP %>% kable(digits = 3,format = "rst") %>% write(file = paste(OutKey,"_reads_distribution.txt",sep=""))


MM %>% fwrite(file = paste(OutKey,"_UpTE_VS_Total_Table_Boxplot.csv",sep = ""))

GENE_ANN %>% fwrite(file = paste(OutKey,"_GENE_Expression_DEG_Annotation.csv",sep = ""))
GENE_ANN[direction=="UP"]$ID %>% write(file = paste(OutKey,"_UP_GENE_Expression_DEG_Annotation.txt",sep = ""))
GENE_ANN[direction=="DOWN"]$ID %>% write(file = paste(OutKey,"_DOWN_GENE_Expression_DEG_Annotation.txt",sep = ""))

TE_RNA %>% fwrite(file = paste(OutKey,"_TE_Expression_DEG_Annotation.csv",sep = ""))
TE_RNA[direction=="UP"] %>% fwrite(file = paste(OutKey,"_Activated_TE_Expression_DEG_Annotation.csv",sep = ""))
TE_RNA[direction=="UP"]$ID %>% write(file = paste(OutKey,"_Activated_TE_ID.txt",sep = ""))



ScatterPlot_TE_Pair %>%  ggsave(filename = paste(OutKey,"_Pair_TE_CPM_Scatterplot.pdf",sep=""),height = 6,width = 6)
GENE_Scatterplot %>% ggsave(filename = paste(OutKey,"_Pair_GENE_RPKM_Scatterplot.pdf",sep=""),height = 6,width = 6)


## Chromosome wide distribution of activated TEs ----
library(gtrellis)

UP<-TE_RNA[direction=="UP"]
UP$log10<-log10(UP$Mutant.mean.cpm - UP$Control.mean.cpm)
GRUP<-GRanges(UP)

Ath_Chromatin<-fread("/data1/linhua/QIANLAB/PROJECT/hic/HOMER/BEDG/Ath_cytoband_Centromere-200K_gieStain.csv")
Ath_Chromatin$group<-0
Ath_Chromatin[gieStain=="acen"]$group<-1
Ath_Chromatin[gieStain=="gneg"]$group<-2
Ath_Chromatin[gieStain=="gpos100"]$group<-3
Ath_Chromatin[gieStain=="gpos95"]$group<-4
Ath_Chromatin[gieStain=="stalk"]$group<-5

col_fun = circlize::colorRamp2(seq(1,5),c("red","white","black","grey","red"))

Chr_df = data.frame(c("Chr1","Chr2","Chr3", "Chr4", "Chr5"), c(1,1,1,1,1),c(30427671,19698289, 23459830, 18585056, 26975502))
Chr_df34 <- data.frame(c("Chr3","Chr4"), c(1,1),c(23459830,18585056))

GG<-GRUP[GRUP$Transposon_Super_Family=="LTR/Gypsy"]
CC<-GRUP[GRUP$Transposon_Super_Family=="LTR/Copia"]
OT<-GRUP[!(GRUP$Transposon_Super_Family %in% c("LTR/Copia","LTR/Gypsy")) & !(GRUP$Class %in% "DNAtransposon")]
DNA<-GRUP[GRUP$Class=="DNAtransposon"]

## Genome-Wide TE distribution ----

TE_Annotation <- fread("/data1/linhua/QIANLAB/PROJECT/te/TE_Annotation-2018-11-17.csv")
table(TE_Annotation$Transposon_Super_Family)
COPIA<-GRanges(TE_Annotation[Transposon_Super_Family=="LTR/Copia"])

Copia_density<-circlize::genomicDensity(region = TE_Annotation[Transposon_Super_Family=="LTR/Copia"][,c("seqnames","start","end")],window.size = 1e5,overlap = F)

Gypsy_Density<-circlize::genomicDensity(region = TE_Annotation[Transposon_Super_Family=="LTR/Gypsy"][,c("seqnames","start","end")],window.size = 1e5,overlap = F)

DNA_TE_Density<-circlize::genomicDensity(region = TE_Annotation[Class=="DNAtransposon"][,c("seqnames","start","end")],window.size = 1e5,overlap = F)

## Chr3 and Chr 4 ----
pdf(paste(OutKey,"_TE_distr_in_Chr34.pdf",sep = ""),height = 6,width = 12)

gtrellis_layout(data=Chr_df34,
                n_track=4,
                #track_ylab=c("",""),
                asist_ticks=F,
                title='TE distribution in 5 chromosomes',
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

## Genome-wide TE activation ----
pdf(paste(OutKey,"_TE_distr_in_Genome.pdf",sep = ""),height = 6,width = 12)

gtrellis_layout(data=Chr_df,
                n_track=4,
                #track_ylab=c("",""),
                asist_ticks=F,
                title='TE distribution in 5 chromosomes',
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
