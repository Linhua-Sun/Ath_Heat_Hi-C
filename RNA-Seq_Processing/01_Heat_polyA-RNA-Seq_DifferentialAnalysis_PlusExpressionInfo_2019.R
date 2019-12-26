## Library R packages and set up ploting parameters ----
library(pacman)
suppressWarnings(suppressPackageStartupMessages(
    p_load(
        ComplexHeatmap,
        crayon,
        data.table,
        dendsort,
        DESeq2,
        dplyr,
        edgeR,
        ggplot2,
        ggpubr,
        ggrepel,
        ggsci,
        gmodels,
        hashmap,
        knitr,
        readr,
        rtracklayer,
        scatterplot3d,
        stringr,
        tidyverse,
        tximport
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

## Aim: audo DEG analysis for multiple sample (with 2 replicates)
## Import quant value for each sample ----
setwd("/data1/linhua/QIANLAB/PROJECT/te/HeatTEReanalysis2019-06-03/")

## data structure of heat_polyA_fc_count.txt (Pls skip first line).
## 1	Geneid
## 2	Chr
## 3	Start
## 4	End
## 5	Strand
## 6	Length
## 7	Col-Control-polyA-1_sortbyname.bam
## 8	Col-Control-polyA-2_sortbyname.bam
## 9	Col-Heat-polyA-1_sortbyname.bam
## 10	Col-Heat-polyA-2_sortbyname.bam
## 11	Col-Recovery-polyA-1_sortbyname.bam
## 12	Col-Recovery-polyA-2_sortbyname.bam

heat_fc<-fread("heat_polyA_fc_count.txt",select = c(1,seq(7,12)),skip = 1)
str_replace_all(colnames(heat_fc),
                pattern = "Col-|polyA-|_sortbyname.bam",
                replacement = "") %>%
    str_replace_all(pattern = "Geneid", "ID") -> colnames(heat_fc)

## Convert raw rounts into data.frame format (for DESeq2 input) ----
heat_fc_df<-as.data.frame(heat_fc)
rownames(heat_fc_df)<-heat_fc_df$ID
heat_fc_df$ID<-NULL

## featurecounts reads distribution analysis (How many reads have been mapped to TE) ----
Ath_GENE_TE_Clean <- fread(file = "/data1/linhua/QIANLAB/PROJECT/te/Ath_GENE_TE_Clean_Pseudo.csv")

full_join(heat_fc,Ath_GENE_TE_Clean,by="ID") %>%
    as.data.table() -> ANN

melt(ANN,id.vars = c("ID","group","locus_type")) -> MANN

MANN %>% 
    group_by(group,variable) %>%
    summarise(NO_mapped_reads=sum(value)) %>%
    group_by(variable) %>%
    mutate(fraction=NO_mapped_reads/sum(NO_mapped_reads)*100) -> PP

PP %>% kable(digits = 3,format = "rst")

pdf(file = "HeatPolyA_RNA_Seq_ReadsDistribution.pdf",height = 8,width = 6)
PP %>% ggbarplot(x = "variable",
                 y = "fraction",
                 group = "group",
                 fill = "group") %>%
    ggpar(
        ylim = c(0, 1),
        x.text.angle = 45,
        xlab = "",
        ylab = "Fraction of reads distribution (%)"
    ) + border() +
    theme(aspect.ratio = 1.5)
dev.off()

## Import into DEseq2 env ----
## Generate experimental design
coldata<-data.frame(treat=c("Control","Control","Heat","Heat","Recovery","Recovery"))
rownames(coldata)<-colnames(heat_fc_df)

F1<-heat_fc_df
F1[rowSums(F1)>=20,] -> F1 ## Filter by 20

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

## (independent PCA part) PCA analysis of top varied 1000 genes ----
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
pc[pc$group== "Heat",]$color<-pal_aaas(alpha = 1)(9)[2]
pc[pc$group== "Recovery",]$color<-pal_aaas(alpha = 1)(9)[3]

## (independent PCA part) 3D PCA PLOT ----
pdf(file = "HeatPolyA_RNA_Seq_PCA_3D_genes_DEseq2.pdf",height = 7,width = 8)
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

## (independent PCA part) 2D PCA plot (ggplot2) ----
pdf(file = "HeatPolyA_RNA_Seq_PCA_2D_genes_DEseq2.pdf",height = 7,width = 8)
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

H_VS_C<-as.data.table(as.data.frame(results(Unit_dds,c("treat","Heat","Control"))),keep.rownames = "ID")[,"direction":="NoChange"][,"Condition":="Heat_vs_Control"]
R_VS_H<-as.data.table(as.data.frame(results(Unit_dds,c("treat","Recovery","Heat"))),keep.rownames = "ID")[,"direction":="NoChange"][,"Condition":="Recovery_vs_Heat"]

temp<-rbind(H_VS_C,R_VS_H)
temp[padj < 0.05 & log2FoldChange>1]$direction<-"UP"
temp[padj < 0.05 & log2FoldChange<(-1)]$direction<-"DOWN"

fwrite(temp,file = "HeatPolyA_RNA_Seq-DEseq2-rbind.csv")

## NO of UP and DOWN GENE and TE----
temp[str_detect(ID,"TE")] %>% group_by(Condition,direction) %>% summarise(count=n()) %>% as.data.table() %$% .[,type:="TE"] -> polyA_DE_TE_No
temp[!str_detect(ID,"TE")] %>% group_by(Condition,direction) %>% summarise(count=n())%>% as.data.table() %$% .[,type:="PCG"] -> polyA_DE_GENE_No
MM<-rbind(polyA_DE_GENE_No,polyA_DE_TE_No)

pdf(file = "polyA-DE-TE-PCG-No-Frac.pdf",height = 8,width = 6)
ggbarplot(data = MM[Condition=="Heat_vs_Control" & direction!="NoChange"],
    x = "type",
    y = "count",
    color = "black",
    fill = "direction",
    palette = pal_aaas()(2),
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
    #x.text.angle = 0,
    ylab = "Fraction and number\n of TE and PCG",
    font.x = c(18),
    font.y = c(18),
    font.tickslab = c(15),
    font.legend = c(15),legend.title = ""
) + 
    #border() +
    #scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))+ # https://stackoverflow.com/questions/13701347/force-the-origin-to-start-at-0
    theme(aspect.ratio = 1/0.618,axis.line.x =element_blank()) # golden ratio portrait
dev.off()


## Prepare CPM table from raw counts ----

ExCPM<-as.data.table(cpm(F1),keep.rownames = "ID")

left_join(ExCPM,heat_fc,by="ID",suffix = c(".cpm", ".rawcount")) %>%
    as.data.table() %$%
    .[,`:=`(Control.mean.cpm=(`Control-1.cpm`+`Control-2.cpm`)/2,
         Heat.mean.cpm=(`Heat-1.cpm`+`Heat-2.cpm`)/2,
         Recovery.mean.cpm=(`Recovery-1.cpm`+`Recovery-2.cpm`)/2
         ),by="ID"][] -> AA

## Join CPM table and DE table for each locus ----
temp[,c("ID","direction","Condition")][,value:="+"][] %$%
    dcast(., ID ~ Condition,value.var="direction") %>% 
    full_join(AA,by = "ID") %>%
    as.data.table() %$% 
    .[,`:=`(Control.mean.cpm.log2=log2(Control.mean.cpm+1),
            Heat.mean.com.log2=log2(Heat.mean.cpm+1),
            Recovery.mean.cpm.log2=log2(Recovery.mean.cpm+1)) ,by="ID"][] -> IIIII

#fwrite(IIIII,file = "Updated_HeatPolyA_RNA_Seq-table.csv")

## Annotation of PCG  ----
#IIIII<-fread("Updated_HeatPolyA_RNA_Seq-table.csv")
GFF<-as.data.table(import("/data1/linhua/QIANLAB/Araport_TOPHAT/Araport11_GFF3_genes_transposons.201606.gff"))
GFF[type=='gene' | type=="transposable_element_gene" |type=="pseudogene"][,c("type","locus_type","ID","Note","symbol","Alias","full_name","seqnames","start","end","width","strand")] -> LLL
LLL$Alias %>% sapply(paste, collapse = "|",simplify = T) -> LLL$Alias 
LLL$Note %>% sapply(paste, collapse = "|",simplify = T) -> LLL$Note


LENGTH <- fread(file = "/data1/linhua/QIANLAB/PROJECT/te/Ath_New_GENE_TE_Length.csv")
left_join(IIIII,LENGTH,by="ID") %>% as.data.table() -> IIIII_ANN

IIIII_ANN[, "Control-1.rpkm":=rpkm(`Control-1.rawcount`,length)]
IIIII_ANN[, "Control-2.rpkm":=rpkm(`Control-2.rawcount`,length)]
IIIII_ANN[, "Heat-1.rpkm":=rpkm(`Heat-1.rawcount`,length)]
IIIII_ANN[, "Heat-2.rpkm":=rpkm(`Heat-2.rawcount`,length)]
IIIII_ANN[, "Recovery-1.rpkm":=rpkm(`Recovery-1.rawcount`,length)]
IIIII_ANN[, "Recovery-2.rpkm":=rpkm(`Recovery-2.rawcount`,length)]

IIIII_ANN[, "Control.mean.rpkm":= (`Control-1.rpkm`+`Control-2.rpkm`)/2,by="ID"]
IIIII_ANN[, "Heat.mean.rpkm":= (`Heat-1.rpkm`+`Heat-2.rpkm`)/2,by="ID"]
IIIII_ANN[, "Recovery.mean.rpkm":= (`Recovery-1.rpkm`+`Recovery-2.rpkm`)/2,by="ID"]

IIIII_ANN[, "Control.mean.rpkm.log2":= log2(Control.mean.rpkm+1),by="ID"]
IIIII_ANN[, "Heat.mean.rpkm.log2":= log2(Heat.mean.rpkm+1),by="ID"]
IIIII_ANN[, "Recovery.mean.rpkm.log2":= log2(Recovery.mean.rpkm+1),by="ID"]

inner_join(IIIII_ANN[!str_detect(ID,"TE")],LLL,by = "ID") %>% as.data.table() -> GENE_ANN

GENE_ANN[order(type)][,c("ID",  "Heat_vs_Control",  "Recovery_vs_Heat",  "Control.mean.cpm",  "Heat.mean.cpm",  "Recovery.mean.cpm",  "type",  "locus_type",  "Note",  "symbol",  "Alias",  "full_name",  "seqnames",  "start",  "end",  "width",  "strand",  "Control.mean.cpm.log2",  "Heat.mean.com.log2",  "Recovery.mean.cpm.log2",  "Control-1.cpm",  "Control-2.cpm",  "Heat-1.cpm",  "Heat-2.cpm",  "Recovery-1.cpm",  "Recovery-2.cpm",  "Control-1.rawcount",  "Control-2.rawcount",  "Heat-1.rawcount",  "Heat-2.rawcount",  "Recovery-1.rawcount",  "Recovery-2.rawcount","Control-1.rpkm","Control-2.rpkm","Heat-1.rpkm","Heat-2.rpkm","Recovery-1.rpkm","Recovery-2.rpkm","Control.mean.rpkm","Heat.mean.rpkm","Recovery.mean.rpkm","Control.mean.rpkm.log2","Heat.mean.rpkm.log2","Recovery.mean.rpkm.log2"),with=F] -> GENE_ANN 

GENE_ANN %>% fwrite("Updated_Heat_Stress_polyA-RNA-Seq_DEG_Ann_Results.csv")

## Annotation of TE  ----
TE_Annotation <- fread("/data1/linhua/QIANLAB/PROJECT/te/TE_Annotation-2019-06-05.csv")

inner_join(IIIII[str_detect(ID, "TE")] , TE_Annotation, by="ID") %>% as.data.table() -> TE_RNA

TE_RNA[,c("ID","Heat_vs_Control","Recovery_vs_Heat","Control.mean.cpm","Heat.mean.cpm","Recovery.mean.cpm","length","Position","seqnames","start","end","Ol_TEGs","Chromatin","Transposon_Family","Transposon_Super_Family","Class","score.GC","score.H3K9me2","Control-1.rawcount","Control-2.rawcount","Heat-1.rawcount","Heat-2.rawcount","Recovery-1.rawcount","Recovery-2.rawcount","Control-1.cpm","Control-2.cpm","Heat-1.cpm","Heat-2.cpm","Recovery-1.cpm","Recovery-2.cpm","Control.mean.cpm.log2","Heat.mean.com.log2","Recovery.mean.cpm.log2")] -> TE_RNA

fwrite(TE_RNA,file = "Updated_Heat_Stress_polyA-RNA-Seq_DEG_TEs_Ann_Results.csv")

#fwrite(TE_RNA[Heat_vs_Control=="UP"],file = "Updated_heat_polyA_UP_TE_Ann.csv")
#fwrite(TE_RNA[Heat_vs_Control=="DOWN"],file = "Updated_heat_polyA_DOWN_TE_Ann.csv")

#TE_RNA[Heat_vs_Control=="UP"][,intersect(colnames(TE_Annotation),colnames(TE_RNA)),with=F][,"group":="Up_TE"]-> TE_ANN
#TE_Annotation[,intersect(colnames(TE_Annotation),colnames(TE_RNA)),with=F][,"group":="Total"] -> Total_ANN

## Heat activated and repressed TEs latest version scatterplot ----
TE_RNA <- fread(file = "Updated_Heat_Stress_polyA-RNA-Seq_DEG_TEs_Ann_Results.csv")

#rbind(TE_RNA[Heat_vs_Control=="UP"][Transposon_Family %in% c("ATCOPIA28","ROMANIAT5","ATCOPIA78")],
#      TE_RNA[Heat_vs_Control=="DOWN"][Transposon_Family %in% c("ATCOPIA93")]) -> Markers

pdf(file = "TE_Heat_VS_Control_polyA_ONSEN_Green_label.pdf",height = 6,width = 6)
ggplot(TE_RNA, aes(x = Control.mean.cpm.log2, y = Heat.mean.com.log2)) +
    geom_point(
        data = TE_RNA[Heat_vs_Control == "NoChange"],
        aes(x = Control.mean.cpm.log2, y = Heat.mean.com.log2),
        color = "black",
        alpha = 0.7,
        size = 3,
        fill = 'lightgrey',
        shape = 21
    ) +
    geom_point(
        data = TE_RNA[Heat_vs_Control == "UP"],
        aes(x = Control.mean.cpm.log2, y = Heat.mean.com.log2),
        color = "black",
        alpha = 0.6,
        size = 3,
        fill = 'red',
        shape = 21
    ) +
    geom_point(
        data = TE_RNA[Heat_vs_Control == "DOWN"],
        aes(x = Control.mean.cpm.log2, y = Heat.mean.com.log2),
        color = "black",
        alpha = 0.6,
        size = 3,
        fill = 'blue',
        shape = 21
    ) +
    geom_point(
        data = TE_RNA[ID %in% c("AT3TE51895","AT1TE43225","AT1TE51360","AT2TE06390","AT4TE10320") ],
        aes(x = Control.mean.cpm.log2, y = Heat.mean.com.log2),
        color = "black",
        alpha = 1,
        size = 4,
        fill = 'green',
        shape = 24
    ) +
    geom_text_repel(
        data = TE_RNA[ID %in% c("AT3TE51895","AT1TE43225","AT1TE51360","AT2TE06390","AT4TE10320")],
        aes(x = Control.mean.cpm.log2, y = Heat.mean.com.log2,label = ID),
        color = "blue"
    )+
    geom_point(
        data = TE_RNA[Transposon_Family %in% c("ATCOPIA78") & length > 4000],
        aes(x = Control.mean.cpm.log2, y = Heat.mean.com.log2),
        color = "black",
        alpha = 1,
        size = 4,
        fill = 'green',
        shape = 24
    ) +
    geom_text_repel(
        data = TE_RNA[Transposon_Family %in% c("ATCOPIA78") & length > 4000],
        aes(x = Control.mean.cpm.log2, y = Heat.mean.com.log2,label = ID),
        color = "blue"
    )+
    scale_x_continuous(limits = c(0, 9),breaks = c(0,2,4,6,8), name = "Control (log2(cpm))") +
    scale_y_continuous(limits = c(0, 9),breaks = c(0,2,4,6,8), name = "Heat (log2(cpm))") +
    geom_abline(intercept = 0, color = "purple",size=1)
dev.off()

