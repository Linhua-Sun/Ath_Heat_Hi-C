## Use HOMER to calculate interaction frequency of KEEs
#```
#(base) [linhua@localhost /data1/linhua/QIANLAB/PROJECT/hic/HOMER/KEE_Pair_Interaction]$cat ../Batch_Cal_KEE_Pairs.sh|grep "KEE7"|grep "KEE8"
#analyzeHiC HeatTagDir/ -res 1000 -window 2000 -ihskb -pos Chr4:11050000-11206537  -pos2 Chr4:15441465-15537500	> KEE7_KEE8_Heat_1K.txt
#analyzeHiC ControlTagDir/ -res 1000 -window 2000 -ihskb -pos Chr4:11050000-11206537  -pos2 Chr4:15441465-15537500	> KEE7_KEE8_Control_1K.txt
#(base) [linhua@localhost /data1/linhua/QIANLAB/PROJECT/hic/HOMER/KEE_Pair_Interaction]$cat ../Batch_Cal_KEE_Pairs.sh|grep "KEE3"|grep "KEE4"
#analyzeHiC HeatTagDir/ -res 1000 -window 2000 -ihskb -pos Chr3:1950000-1971581    -pos2 Chr3:3100000-3121455	> KEE3_KEE4_Heat_1K.txt
#analyzeHiC ControlTagDir/ -res 1000 -window 2000 -ihskb -pos Chr3:1950000-1971581    -pos2 Chr3:3100000-3121455	> KEE3_KEE4_Control_1K.txt
#```

## ----
library(plot3D)
## Persp3D on HiC selected regions:
setwd("/data1/linhua/QIANLAB/PROJECT/hic/HOMER")
Control<-fread("./KEE_Pair_Interaction/KEE7_KEE8_Control_1K.txt")
## Control[,c("HiCMatrix (directory=ControlTagDir/)","Regions"):=NULL]

Control<-as.data.frame(Control)
rownames(Control)<-paste(round(as.numeric(str_remove_all(Control$Regions,"Chr4-"))/1000/1000,digits = 3),"M",sep="")
colnames(Control)<-paste(round(as.numeric(str_remove_all(colnames(Control),"Chr4-"))/1000/1000,digits = 3),"M",sep="")
Control[,c(1,2)]<-NULL
Control <- as.matrix(Control)


rLab<-as.numeric(str_remove_all(rownames(Control),"M"))
cLab<-as.numeric(str_remove_all(colnames(Control),"M"))

pdf(file = "20191129_Control_KEE7_KEE8_Interaction_persp3D.pdf",height = 6,width =8)
persp3D(
    x =rLab ,y = cLab,
    z = Control,
    theta = 50,
    phi = 25,
    shade = 0.1,
    ticktype = "detailed",
    expand = 1,
    bty = "b2",
    zlim = c(0, 600),
    clim = c(0, 500),
    xlab = "",
    ylab = "",
    zlab = "Interaction Score",
    axes = T
)
graphics.off()

## ----
Heat<-fread("./KEE_Pair_Interaction/KEE7_KEE8_Heat_1K.txt")
## HeatTagDir[,c("HiCMatrix (directory=HeatTagDir/)","Regions"):=NULL]

Heat<-as.data.frame(Heat)
rownames(Heat)<-paste(round(as.numeric(str_remove_all(Heat$Regions,"Chr4-"))/1000/1000,digits = 3),"M",sep="")
colnames(Heat)<-paste(round(as.numeric(str_remove_all(colnames(Heat),"Chr4-"))/1000/1000,digits = 3),"M",sep="")
Heat[,c(1,2)]<-NULL
Heat <- as.matrix(Heat)


rLab<-as.numeric(str_remove_all(rownames(Heat),"M"))
cLab<-as.numeric(str_remove_all(colnames(Heat),"M"))

pdf(file = "20191129_Heat_KEE7_KEE8_Interaction_persp3D.pdf",height = 6,width =8)
persp3D(
    x =rLab ,y = cLab,
    z = Heat,
    theta = 50,
    phi = 25,
    shade = 0.1,
    ticktype = "detailed",
    expand = 1,
    bty = "b2",
    zlim = c(0, 600),
    clim = c(0, 500),
    xlab = "",
    ylab = "",
    zlab = "Interaction Score",
    axes = T
)
graphics.off()
