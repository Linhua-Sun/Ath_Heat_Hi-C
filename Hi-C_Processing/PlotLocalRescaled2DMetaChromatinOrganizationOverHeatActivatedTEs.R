## Import R packages and Define R functions ----
library(pacman)
suppressWarnings(suppressPackageStartupMessages(
    p_load(
        data.table,
        ComplexHeatmap
    )
))


setwd("/data1/linhua/QIANLAB/PROJECT/hic/IDE_201911/")

CCal<- function (file,Name) {
    PP<-fread(file)
    summary(unlist(PP))
    Prefix<-str_remove_all(file,".txt")
    print(str_remove_all(file,".txt|1KBal|_001K_TESm"))
    pdf(file = paste(Prefix,".pdf",sep = ""))
    Heatmap(
        PP,
        col=colorRamp2(c(1,1.15,1.3), c("blue", "white", "red")),
        cluster_columns = F,
        cluster_rows = F,
        show_row_names = F,
        show_column_names = F,
        width = unit(3, "cm"),height = unit(3, "cm"),
        column_title_side = "top",
        column_title = Name
    ) -> Temp
    dev.off()
    return(Temp)
    }

## Heat-activated Tes ----

Control_001K<-CCal(file = "1KBalControl_001K_TESm.txt",Name = "Control")
Control_Rep1_001K<-CCal(file = "1KBalControl_Rep1_001K_SmTe.txt",Name = "Control_Rep1")
Control_Rep2_001K<-CCal(file = "1KBalControl_Rep2_001K_SmTe.txt",Name = "Control_Rep2")
Heat_001K<-CCal(file = "1KBalHeat_001K_TESm.txt",Name = "Heat")
Heat_Rep1_001K<-CCal(file = "1KBalHeat_Rep1_001K_SmTe.txt",Name = "Heat_Rep1")
Heat_Rep2_001K<-CCal(file = "1KBalHeat_Rep2_001K_SmTe.txt",Name = "Heat_Rep2")

pdf(file = "HeatActivatedTeMeta2DAPA.pdf",width = 10,height = 8)
Control_001K + 
Control_Rep1_001K + 
Control_Rep2_001K + 
Heat_001K + 
Heat_Rep1_001K + 
Heat_Rep2_001K
dev.off()

