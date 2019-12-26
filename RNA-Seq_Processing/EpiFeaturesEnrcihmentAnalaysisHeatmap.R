## Library R packages ----
library(pacman)
suppressWarnings(suppressPackageStartupMessages(
    p_load(
        ComplexHeatmap,
        stringr,
        data.table,
        dplyr
        
    )
))
library(circlize)
## Starting ----

Args<-commandArgs(trailingOnly = T)

if (length(Args)!=1) {
    stop(
        "At least 1 argument must be supplied: Rscriipt script.R  test_path.txt",
        call. = FALSE
    )
}

Args[1]<-"HHH.txt"

Inp_Files<-readLines(Args[1])

T1<-fread(Inp_Files[1],header=T)

T1[,"Class":=unlist(str_split(Inp_Files[1],pattern = "-"))[1]]

for (i in 2: length(Inp_Files)) {
    TEMP<-fread(Inp_Files[i],header=T)
    TEMP[,"Class":=unlist(str_split(Inp_Files[i],pattern = "-"))[1]]
    T1<-rbind(T1,TEMP)
}

LEVELS<-unique(T1$Class)

## STEP1 ----
TT1<-T1[Class==LEVELS[1]]
compare_means(formula = score ~ group , data = TT1,ref.group = "Total") %>% as.data.table() -> LL

LL$GROUP<-LEVELS[1]
LL$DIFF<-0
unique(TT1$group)-> GG
for (i in 2:length(GG)) {
    log2(mean(na.omit(TT1[group==GG[i]]$score))/ mean(na.omit(TT1[group=="Total"]$score)) ) -> LL[group2==GG[i]]$DIFF
}

## STEP2 ----

for (i in 2:length(LEVELS)) {
    TEMP<-T1[Class==LEVELS[i]]
    compare_means(formula = score ~ group , data = TEMP,ref.group = "Total") %>% as.data.table() -> TT
    TT$GROUP<-LEVELS[i]
    TT$DIFF<-0
    unique(TEMP$group)-> GG
    for (j in 2:length(GG)) {
        log2(mean(na.omit(TEMP[group==GG[j]]$score))/ mean(na.omit(TEMP[group=="Total"]$score)) ) -> TT[group2==GG[j]]$DIFF
    }
    LL<-rbind(LL,TT)
}

## Generate Heatmap ----
LL[,c("group2","p.adj","GROUP","DIFF")] -> QQ

QQ -> PP

PP[p.adj>0.05]$DIFF<-NA

PP %>% dcast(formula = GROUP ~ group2) -> HH

HH<-as.data.frame(HH)
rownames(HH)<-HH$GROUP
HH$GROUP<-NULL


NewOrder<-c("H2A.W","H2A.X","H2A.Z","H3.1","H3.3","H3K9ac","H4K16ac","H3K14ac","H3K23ac","H3K27ac","H3K36ac","H3K56ac","H3K36me3","H3K4me1","H3K4me2","H3K4me3","H3K9me2","H3K27me1","H3K27me3")

pdf(
    file = str_c("TEMP",str_remove_all(Args[1], ".txt"), "_Heatmap_Histones.pdf", sep = ""),
    height = 8,
    width = 6
)
HH  %>%  Heatmap(
    cluster_rows = F,col = colorRamp2(c(-0.5, 0, 2), c("blue", "white", "red")),
    cluster_columns = F,
    row_order = NewOrder,
    rect_gp = gpar(col = "white",
                   lty = 1,
                   lwd = 0.5),
    na_col = "darkgrey",name = "Score",
    column_names_rot = 30,
    border = F,
    show_row_names = T,
    width = unit(3, "cm"),
    height = unit(9, "cm")
)
dev.off()