#!/usr/bin/env Rscript
## Library R packages ----
library(pacman)
suppressWarnings(suppressPackageStartupMessages(
    p_load(
        crayon,
        data.table,
        dplyr,
        ggplot2,
        ggpubr,
        ggsci,
        stringr,
        tximport,
        rtracklayer,
        ggrepel,
        viridis
    )
))

source("/data1/linhua/QIANLAB/R_Functions/theme_linhua.R")

theme_set(theme_linhua(base_size = 18,base_family = "",legend = "right"))
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
    aspect.ratio = 1
)
## ----

#setwd("/data1/linhua/QIANLAB/PROJECT/hic/HOMER")

Args<-commandArgs(trailingOnly = T)

if (length(Args)==0) {
    stop("At least 1 argument must be supplied: Rscriipt script.R TEST_TE.csv", call.=FALSE)
}

#for test:
#Args<-NULL
#Args[1]<-c("/data1/linhua/QIANLAB/PROJECT/hic/HOMER/Control_PC1_R20KW60K_ByArmsH3K4.PC1.bedGraph")
#Args[2]<-c("/data1/linhua/QIANLAB/PROJECT/hic/HOMER/Heat_PC1_R20KW60K_ByArmsH3K4.PC1.bedGraph")
#BEDG1<-as.data.table(import(Args[1]))[,"coordinate":=paste(seqnames,":",start,end,sep = "")]
#BEDG2<-as.data.table(import(Args[2]))[,"coordinate":=paste(seqnames,":",start,end,sep = "")]

ImportHomerBedGraph <- function(FileInput) {
    data1<-fread(FileInput,skip = 1)
    colnames(data1)<-c("seqnames","start","end","score")
    return(data1[,"coordinate":=paste(seqnames,":",start,"-",end,sep = "")][])
}

BEDG1<-ImportHomerBedGraph(FileInput = Args[1])
BEDG2<-ImportHomerBedGraph(FileInput = Args[2])

FileNamePrefix<-str_c(str_remove_all(basename(Args[1]),".bedGraph|_Res20KWin60K_SplitByCentromereH3K4me3"),"-VS-",str_remove_all(basename(Args[2]),".bedGraph|_Res20KWin60K_SplitByCentromereH3K4me3"),"_PC1",sep="")

FileNamePrefix <- "CheckABSwitch"

#R10KW10K<-as.data.table(import("Control_PC1_R10KW10K.PC1.bedGraph"))[,"coordinate":=paste(seqnames,":",start,end,sep = "")]
#R10KW20K<-as.data.table(import("Control_PC1_R10KW20K.PC1.bedGraph"))[,"coordinate":=paste(seqnames,":",start,end,sep = "")]

inner_join(BEDG1,BEDG2,by="coordinate") %>% as.data.table() -> MM

RES<-cor.test(MM$score.x,MM$score.y)

write(paste("Cor: ",RES$estimate),file = paste(FileNamePrefix,"_Cor.txt",sep=""))

Chr1_Cor<-paste("Chr1_Cor:",cor.test(MM[seqnames.x=="Chr1"]$score.x,MM[seqnames.x=="Chr1"]$score.y)$estimate)
Chr2_Cor<-paste("Chr2_Cor:",cor.test(MM[seqnames.x=="Chr2"]$score.x,MM[seqnames.x=="Chr2"]$score.y)$estimate)
Chr3_Cor<-paste("Chr3_Cor:",cor.test(MM[seqnames.x=="Chr3"]$score.x,MM[seqnames.x=="Chr3"]$score.y)$estimate)
Chr4_Cor<-paste("Chr4_Cor:",cor.test(MM[seqnames.x=="Chr4"]$score.x,MM[seqnames.x=="Chr4"]$score.y)$estimate)
Chr5_Cor<-paste("Chr5_Cor:",cor.test(MM[seqnames.x=="Chr5"]$score.x,MM[seqnames.x=="Chr5"]$score.y)$estimate)
c(Chr1_Cor,
  Chr2_Cor,
  Chr3_Cor,
  Chr4_Cor,
  Chr5_Cor) %>% write(file = paste(FileNamePrefix,"_CorByChr.txt",sep=""))
## ----

AB<-nrow(MM[score.x>0 & score.y<0]) ## A -> B AB
BA<-nrow(MM[score.x<0 & score.y>0]) ## B -> A BA
AA<-nrow(MM[score.x>=0 & score.y>=0]) ## A -> A AA
BB<-nrow(MM[score.x<=0 & score.y<=0]) ## B -> B BB

data.table(N = c(AB,
                 BA,
                 AA,
                 BB),
           group = c("A -> B",
                     "B -> A",
                     "A -> A",
                     "B -> B")) -> df

df %>% mutate(lab= paste(round(N/sum(N)*100,1),"%",sep = "")) -> df

df$group <- factor(c("A -> B",
                     "B -> A",
                     "A -> A",
                     "B -> B"),
                   levels = c("A -> B",
                              "B -> A",
                              "A -> A",
                              "B -> B"))

pdf(paste(FileNamePrefix,"AB_PieChart.pdf",sep = ""))
ggpie(df, "N", label = "lab",lab.font = "white",
      fill = "group", color = "white",palette = "aaas")
dev.off()
dt<-as.data.table(df)

## ----
XName<-str_c(str_remove_all(basename(Args[1]),".bedGraph|_PC1_R20KW60K_ByArmsH3K4.PC1")," PC1",sep = "")
YName<-str_c(str_remove_all(basename(Args[2]),".bedGraph|_PC1_R20KW60K_ByArmsH3K4.PC1")," PC1",sep = "")

annotations <- data.frame(
    xpos =  c(-Inf,Inf),
    ypos =  c( Inf,-Inf),
    annotateText = c(
        paste("B -> A : ", dt[group == "B -> A"]$N, " (", dt[group == "B -> A"]$lab, ")", sep = ""),
        paste("A -> B : ", dt[group == "A -> B"]$N, " (", dt[group == "A -> B"]$lab, ")", sep = "")
    ), 
    hjustvar = c(-0.5,1) ,
    vjustvar = c(1,-0.5))

pdf(file = paste(FileNamePrefix,"_Scatterplot.pdf",sep = ""),width = 8,height = 6)

P <-ggplot(MM, aes(x = score.x, y = score.y)) +
    scale_x_continuous(name = XName) +
    scale_y_continuous(name = YName)

L1 <- P + 
    geom_point(color="lightgrey",shape=21,size=1)+ 
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText,color=annotateText),color=c("red","blue"))

L2 <- geom_point(data = MM[score.x < 0 & score.y > 0 ],aes(x = score.x, y = score.y),color = "red",alpha = 1,size = 2,shape=21)
L3 <- geom_point(data = MM[score.x > 0 & score.y < 0 ],aes(x = score.x, y = score.y),color = "blue",alpha = 1,size = 2,shape=21)

L1 + L2 + L3 + 
    stat_cor(method = "pearson") +    
    geom_vline(xintercept = 0,col="black",size=0.6,linetype="solid")+    
    geom_hline(yintercept = 0,col="black",size=0.6,linetype="solid")+    
    geom_abline(intercept = 0,col="black",size=0.6,linetype="solid",slope = 1)
dev.off()


## ----

#GRanges(BEDG1[,c(1:5)]) %>% export.bed(con = "All_PCA_UsedBinsPer20K.bed")

MM[score.x>0 & score.y<0][,c(1:5)] -> MM_AB
colnames(MM_AB)<-str_remove_all(string = colnames(MM_AB),pattern = ".x")
MM_AB$seqnames<-str_replace_all(MM_AB$seqnames,pattern = "Chr",replacement = "chr")
GRanges(MM_AB) %>% export.bed(con = paste(FileNamePrefix,"_A-To-B.bed",sep = ""))

MM[score.x<0 & score.y>0]-> MM_BA
colnames(MM_BA)<-str_remove_all(string = colnames(MM_BA),pattern = ".x")
MM_BA$seqnames<-str_replace_all(MM_BA$seqnames,pattern = "Chr",replacement = "chr")
GRanges(MM_BA) %>% export.bed(con = paste(FileNamePrefix,"_B-To-A.bed",sep = ""))

