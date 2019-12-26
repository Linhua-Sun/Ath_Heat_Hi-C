## Set Up rowScale function ----
rowScale = function(x,
                    center = TRUE,
                    scale = TRUE,
                    add_attr = TRUE,
                    rows = NULL,
                    cols = NULL) {
    library(matrixStats)
    
    if (!is.null(rows) && !is.null(cols)) {
        x <- x[rows, cols, drop = FALSE]
    } else if (!is.null(rows)) {
        x <- x[rows, , drop = FALSE]
    } else if (!is.null(cols)) {
        x <- x[, cols, drop = FALSE]
    }
    
    ################
    # Get the column means
    ################
    cm = rowMeans(x, na.rm = TRUE)
    ################
    # Get the column sd
    ################
    if (scale) {
        csd = rowSds(x, center = cm)
    } else {
        # just divide by 1 if not
        csd = rep(1, length = length(cm))
    }
    if (!center) {
        # just subtract 0
        cm = rep(0, length = length(cm))
    }
    x = (x - cm) / csd
    if (add_attr) {
        if (center) {
            attr(x, "scaled:center") <- cm
        }
        if (scale) {
            attr(x, "scaled:scale") <- csd
        }
    }
    return(x)
}

## Import TE expression ----
TE_RNA <- fread(file = "Updated_Heat_Stress_polyA-RNA-Seq_DEG_TEs_Ann_Results.csv")
colnames(TE_RNA)<-str_replace_all(colnames(TE_RNA),"com","cpm")
TE_RNA[,Path:=paste(Heat_vs_Control,Recovery_vs_Heat,sep = ":")]
TE_RNA[,c("ID",
          "Heat_vs_Control",
          "Recovery_vs_Heat",
          "Control.mean.cpm",
          "Heat.mean.cpm",
          "Recovery.mean.cpm",
          "seqnames",
          "length",
          "Chromatin",
          "Transposon_Family",
          "Transposon_Super_Family",
          "Class",
          "score.GC",
          "score.H3K9me2",
          "Control.mean.cpm.log2",
          "Heat.mean.cpm.log2",
          "Recovery.mean.cpm.log2",
          "Path"),with=F] -> TE
TE$average.cpm.log2 <- rowMeans(TE[,c("Control.mean.cpm.log2","Heat.mean.cpm.log2","Recovery.mean.cpm.log2")])

TE[Path!="NoChange:NoChange"][order(Path,decreasing = T)] -> Full_Based
Full_Based[,group:=as.numeric(factor(Path))]

temp<-rowScale(Full_Based[,c("Control.mean.cpm.log2","Heat.mean.cpm.log2","Recovery.mean.cpm.log2")]) # ## https://www.r-bloggers.com/a-faster-scale-function/
colnames(temp)<-paste(colnames(temp),".z_score",sep="")
Full_Based_Extend<-cbind(Full_Based,temp)

# Convert to data.frame ----
setnames(Full_Based_Extend,c("Control.mean.cpm.log2.z_score","Heat.mean.cpm.log2.z_score","Recovery.mean.cpm.log2.z_score"),c("Control","Heat","Recovery"))

Full_Based_Extend_DF<-data.frame(Full_Based_Extend)
rownames(Full_Based_Extend_DF)<-Full_Based_Extend_DF$ID
Full_Based_Extend_DF$ID<-NULL

## Set up Chromatin_Col_Key ----
# c("CA"="gainsboro",
# "PR"="dimgrey",
# "knob"="blue",
# "centromere"="black",
# "KEE"="red") -> context_Col_Key

 c("8"=pal_aaas(alpha = 0.7)(8)[8],
 "7"=pal_aaas(alpha = 0.7)(8)[7],
 "6"=pal_aaas(alpha = 0.7)(8)[6],
 "5"=pal_aaas(alpha = 0.7)(8)[5],
 "4"=pal_aaas(alpha = 0.7)(8)[4],
 "3"=pal_aaas(alpha = 0.7)(8)[3],
 "2"=pal_aaas(alpha = 0.7)(8)[2],
 "1"=pal_aaas(alpha = 0.7)(8)[1]) -> group_Col_Key

c("Chr1"="black",
  "Chr2"="white",
  "Chr3"="black",
  "Chr4"="white",
  "Chr5"="black") -> group_seqnames_Key

## Heatmap Annotation preparation ----
Annotation_Rows <-
    rowAnnotation(
        df = Full_Based_Extend_DF[, c("average.cpm.log2", "Class", "seqnames")],
        col = list(
            #group = group_Col_Key,
            Class = c(
                "Unknown" = "grey",
                "retrotransposon" = "#EE00007F",
                "DNAtransposon" = "lightblue"
            ),
            #context = context_Col_Key,
            seqnames=group_seqnames_Key,
            average.cpm.log2 = colorRamp2(c(0, 6), c("white", "red"))
        )
    )

Clean_Annotation_Rows <-
    rowAnnotation(
        df = Full_Based_Extend_DF[, c("average.cpm.log2","group")],
        col = list(group = group_Col_Key,
                   average.cpm.log2 = colorRamp2(c(0, 6), c("white", "red")))
        
    )

 ## AllChangedTE-Row-Hclust-Z-score-Heatmap ----

pdf(file = "AllChangedTE_ChrSortedZ-scoreHeatmap.pdf",height = 10,width = 6)
Full_Based_Extend_DF[, (c("Control", "Heat", "Recovery"))] %>%
    Heatmap(
        name = "Z-score",
        #col = colorRamp2(c(-1, 0, 1), c("green", "white", "red")),
        use_raster = F,
        cluster_columns = F,
        cluster_rows = T,
        #clustering_distance_rows = "pearson",
        split = as.numeric(factor(Full_Based_Extend_DF$Path)),
        column_title = "TE",
        show_row_dend = T,
        show_row_names = F,
        rect_gp = gpar(col = NA, lty = 2, lwd = 0,border=F),
        column_title_gp = gpar(fontsize = 10),
        width = unit(5, "cm")
        #rect_gp = gpar(col = "white",lwd=0.25),
    )  + Clean_Annotation_Rows
dev.off()

## Only Up TE (Heat VS Control) ----

DF<-as.data.frame(Full_Based_Extend[group %in% c(6,7)] )

pdf(file = "HeatActivatedTE_ChrSortedZ-scoreHeatmap.pdf",height = 12,width = 8)
DF[, rev(c("Control", "Heat", "Recovery"))] %>%
    Heatmap(
        name = "Z-score",
        #col = colorRamp2(c(-1, 0, 1), c("green", "white", "red")),
        cluster_columns = F,
        cluster_rows = F,
        #clustering_distance_rows = "pearson",
        split = as.numeric(factor(DF$Path)),
        column_title = "TE",
        show_row_dend = F,
        show_row_names = F,
        column_title_gp = gpar(fontsize = 10),
        width = unit(3.5, "cm"),use_raster = F
        #rect_gp = gpar(col = "white",lwd=0.25),
    ) + Annotation_Rows
dev.off()

## Not Now ----
## Boxplot analysis (per group) ----
## pdf(file = "Heat_Activated_TEs_Group1.pdf",height = 8,width = 6)
## Full_Based_Extend_DF[Full_Based_Extend_DF$group==6,][,c("Control.mean.cpm.log2","Heat.mean.cpm.log2","Recovery.mean.cpm.log2")] %>% melt() -> TTT
## 
## Stats<-compare_means(formula = value ~ variable , data = TTT)
## ggplot(TTT, aes(x =variable , y = value , fill = variable)) +
##     stat_boxplot(geom = 'errorbar', width = 0.5) +
##     geom_boxplot(outlier.colour = "black",width = 0.5,outlier.size = 1,notch = F,notchwidth = 0,outlier.shape = 1) +
##     stat_summary(fun.y = median,geom = "line",aes(group = 1),size = 1,color = "deepskyblue2") +
##     stat_summary(fun.y = median,geom = "point",shape = 21,size = 3,fill = "yellow",alpha = 0.8) +
##     scale_fill_manual(values = pal_aaas()(9)[3:1],guide=FALSE) +
##     #labs(x = paste(Stats$method," p-value: ",Stats$p.format,sep = ""), y = "Value") +
##     labs(x = "",y = "Value (log2)") +
##     scale_y_continuous(expand = c(0, 0),breaks = scales::pretty_breaks(n = 5))+
##     scale_x_discrete(labels = c("Control","Heat","Recovery"))+
##     border() +theme(aspect.ratio = 2)
## #stat_compare_means(aes(label = paste0("p ", ..p.format..)),size=5)
## dev.off()
## 
## 
## pdf(file = "Heat_Activated_TEs_Group2.pdf",height = 8,width = 6)
## Full_Based_Extend_DF[Full_Based_Extend_DF$group==7,][,c("Control.mean.cpm.log2","Heat.mean.cpm.log2","Recovery.mean.cpm.log2")] %>% melt() -> TTT
## 
## Stats<-compare_means(formula = value ~ variable , data = TTT)
## ggplot(TTT, aes(x =variable , y = value , fill = variable)) +
##     stat_boxplot(geom = 'errorbar', width = 0.5) +
##     geom_boxplot(outlier.colour = "black",width = 0.5,outlier.size = 1,notch = F,notchwidth = 0,outlier.shape = 1) +
##     stat_summary(fun.y = median,geom = "line",aes(group = 1),size = 1,color = "deepskyblue2") +
##     stat_summary(fun.y = median,geom = "point",shape = 21,size = 3,fill = "yellow",alpha = 0.8) +
##     scale_fill_manual(values = pal_aaas()(9)[3:1],guide=FALSE) +
##     #labs(x = paste(Stats$method," p-value: ",Stats$p.format,sep = ""), y = "Value") +
##     labs(x = "",y = "Value (log2)") +
##     scale_y_continuous(expand = c(0, 0),breaks = scales::pretty_breaks(n = 5))+
##     scale_x_discrete(labels = c("Control","Heat","Recovery"))+
##     border() +theme(aspect.ratio = 2)
## #stat_compare_means(aes(label = paste0("p ", ..p.format..)),size=5)
## dev.off()