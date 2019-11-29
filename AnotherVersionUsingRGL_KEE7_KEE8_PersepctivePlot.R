library(rgl)
library(plot3D)
library(ggsci)
library(scales)
library(pheatmap)
library(plot3D)
library(stringr)

## https://2heng.xin/2017/11/19/create-3d-plot-colored-according-to-the-z-axis/
## https://stackoverflow.com/questions/40577408/3d-plot-in-r-using-persp3d-axis-issues
library(data.table)

LL<-fread("~/Downloads/K7K8SM_J.txt")

LL[,i:=paste(LL$V1,"_",LL$V2,"_",LL$V3,sep = "")]
LL[,j:=paste(LL$V4,"_",LL$V5,"_",LL$V6,sep = "")]

A<-1:length(unique(LL$i))
B<-1:length(unique(LL$j))

m0<-matrix(data = 0,nrow = length(unique(LL$i)) , ncol=length(unique(LL$j)))

rownames(m0)<-unique(LL$i)
colnames(m0)<-unique(LL$j)

## View(m0[1:10,1:10])

for (Line in 1:nrow(LL)) {
    print(Line)
    
    print(m0[LL[Line]$i,LL[Line]$j])
    
    m0[LL[Line]$i,LL[Line]$j]<-LL[Line]$V8

}

pheatmap(mat = log2(m0+1),cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)

A<- as.numeric(str_split_fixed(rownames(m0),pattern = "_",n = 3)[,2])
B<- as.numeric(str_split_fixed(colnames(m0),pattern = "_",n = 3)[,2])

BACK<-m0

clab<-sort(as.integer(str_split_fixed(colnames(m0),"_",3)[,2]))
rlab<-sort(as.integer(str_split_fixed(rownames(m0),"_",3)[,2]))
colnames(m0)<-clab
rownames(m0)<-rlab

vertcol <- cut(m0, 100)

persp3d(
    clab,
    rlab,
    m0,
    color =jet.col (n = 100, alpha = 1) [vertcol],
    alpha = 0.7,
    aspect = c(100, 200, 20),
    xlab = "clab",zlim = c(0,10),
    ylab = "rlab",
    zlab = "m0"
)

## Usage ----

jet.col (n = 100, alpha = 1) %>% show_col()
jet2.col (n = 100, alpha = 1) %>% show_col()
gg.col (n = 100, alpha = 1) %>% show_col()
gg2.col (n = 100, alpha = 1) %>% show_col()
ramp.col (col = c("grey", "black"), n = 100, alpha = 1) %>% show_col()
alpha.col (col = "grey", alpha = 0.5) %>% show_col()