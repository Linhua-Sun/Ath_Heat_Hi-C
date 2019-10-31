#!/usr/bin/env Rscript

## Aim: Compare two Hi-C matrix based on HiCPro output
## Author: Linhua Sun Linhua_Sun@pku.edu.cn
## Time: 2019-10-30

## SetUp ----
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(crayon))

option_list = list(
    make_option(
        c("-c", "--Control"),
        type = "character",
        default = NULL,
        help = "dataset file name",
        metavar = "character"
    ),
    make_option(
        c("-t", "--Treatment"),
        type = "character",
        default = NULL,
        help = "dataset file name",
        metavar = "character"
    ),
    make_option(
        c("-b", "--bin"),
        type = "integer",
        default = NULL,
        help = "bin size",
        metavar = "number"
    ),
    make_option(
        c("-C", "--CHR"),
        type = "logical",
        default = T,
        help = "Plot for each chromosome [default= %default]",
        metavar = "logical"
    ),
    make_option(
        c("-o", "--out"),
        type = "character",
        default = "temp",
        help = "output file name prefix [default= %default]",
        metavar = "character"
    )
)

opt_parser = OptionParser(option_list = option_list)

opt = parse_args(opt_parser)

if (is.null(opt$Control)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

cat(yellow(paste("The control file is ", opt$Control,  sep = "")),sep = "\n")
cat(yellow(paste("The treatment file is ", opt$Treatment,  sep = "")),sep = "\n")
cat(yellow(paste("The bin size is ", opt$bin,  sep = "")),sep = "\n")
cat(yellow(paste("The output file prefix is ", opt$out, sep = "")),sep = "\n")

## Library R packages ----
library(pacman)
suppressWarnings(suppressPackageStartupMessages(
    p_load(
        HiCdatR,
        fields,
        fst
    )
))

#source("/data1/linhua/QIANLAB/PROJECT/hic/FPC/BIN100K_Col_Matrix/BU/Modified_Plot_Heatmap_Base.R")
source("/data1/linhua/QIANLAB/PROJECT/hic/Hi-C-2019/HiC_Heatmap_Generator.R")

## Input and Output ----
pathToScripts <- "/data1/linhua/QIANLAB/PROJECT/hic/MWSchmid_test/Rscripts/"
f.source.organism.specific.code(file.path(pathToScripts, "HiCdat-A-thaliana-TAIR10.R"))

## Example ----

#setwd("/data1/linhua/QIANLAB/PROJECT/hic/Hi-C-2019")
#opt<-NULL
#opt$Control<-"Col22_total_sample_100000.matrix"
#opt$Treatment<-"IAA_0h_100000.matrix"
#opt$bin<-100000
#opt$out<-"RootShoot"
#opt$CHR<-T

binSize <- opt$bin

sampleList <- list(
    Control = opt$Control,
    Treatment = opt$Treatment
)

# load all samples into a list, access an individual HiC matrix with binMatList[[sampleName]] or binMatList$sampleName

binMatList <- f.load.samples(getwd(), sampleList, binSize, 50) ## 5959*5959 MATRIX

write.fst(as.data.frame(binMatList$Control),  path = paste(opt$out,"_",opt$Control,"_binMat.fst",sep = ""),compress = 100)
write.fst(as.data.frame(binMatList$Treatment),path = paste(opt$out,"_",opt$Treatment,"_binMat.fst",sep = ""),compress = 100)

## Control and Treatment Raw Normalized Heatmap ----

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
pdf(file = paste(opt$out,"_",str_remove_all(opt$Control,".matrix"),"_WGRawMap.pdf",sep=""),width = 19,height = 19)
Genome.f.plot.XY.matrix(
    matrixToPlot = binMatList$Control,
    binSize = opt$bin,
    axStep = 10000000 ,
    chromA = "ALL",
    startA = 0,
    endA = 0,
    chromB = "ALL",
    startB = 0,
    endB = 0,drawGrid = T
)-> RawIntControl
dev.off()

if (opt$CHR) {
    pdf(file = paste(opt$out, "_", str_remove_all(opt$Control,".matrix"), "_PerChrRawMap.pdf", sep = ""),width = 19,height = 19)
    Chrf.plot.XY.matrix(matrixToPlot = binMatList$Control,binSize = opt$bin,axStep = 5000000 ,chromA = "Chr1",startA = 0,endA = 0,chromB = "Chr1",startB = 0,endB = 0)
    Chrf.plot.XY.matrix(matrixToPlot = binMatList$Control,binSize = opt$bin,axStep = 5000000 ,chromA = "Chr2",startA = 0,endA = 0,chromB = "Chr2",startB = 0,endB = 0)
    Chrf.plot.XY.matrix(matrixToPlot = binMatList$Control,binSize = opt$bin,axStep = 5000000 ,chromA = "Chr3",startA = 0,endA = 0,chromB = "Chr3",startB = 0,endB = 0)
    Chrf.plot.XY.matrix(matrixToPlot = binMatList$Control,binSize = opt$bin,axStep = 5000000 ,chromA = "Chr4",startA = 0,endA = 0,chromB = "Chr4",startB = 0,endB = 0)
    Chrf.plot.XY.matrix(matrixToPlot = binMatList$Control,binSize = opt$bin,axStep = 5000000 ,chromA = "Chr5",startA = 0,endA = 0,chromB = "Chr5",startB = 0,endB = 0)
    dev.off()
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

pdf(file = paste(opt$out,"_",str_remove_all(opt$Treatment,".matrix"),"_WGRawMap.pdf",sep=""),width = 19,height = 19)
Genome.f.plot.XY.matrix(
    matrixToPlot = binMatList$Treatment,
    binSize = opt$bin,
    axStep = 10000000 ,
    chromA = "ALL",
    startA = 0,
    endA = 0,
    chromB = "ALL",
    startB = 0,
    endB = 0,drawGrid = T
)-> RawIntTreatment
dev.off()

if (opt$CHR) {
    pdf(file = paste(opt$out, "_", opt$Treatment, "_PerChrRawMap.pdf", sep = ""),width = 19,height = 19)
    Chrf.plot.XY.matrix(matrixToPlot = binMatList$Treatment,binSize = opt$bin,axStep = 5000000 ,chromA = "Chr1",startA = 0,endA = 0,chromB = "Chr1",startB = 0,endB = 0)
    Chrf.plot.XY.matrix(matrixToPlot = binMatList$Treatment,binSize = opt$bin,axStep = 5000000 ,chromA = "Chr2",startA = 0,endA = 0,chromB = "Chr2",startB = 0,endB = 0)
    Chrf.plot.XY.matrix(matrixToPlot = binMatList$Treatment,binSize = opt$bin,axStep = 5000000 ,chromA = "Chr3",startA = 0,endA = 0,chromB = "Chr3",startB = 0,endB = 0)
    Chrf.plot.XY.matrix(matrixToPlot = binMatList$Treatment,binSize = opt$bin,axStep = 5000000 ,chromA = "Chr4",startA = 0,endA = 0,chromB = "Chr4",startB = 0,endB = 0)
    Chrf.plot.XY.matrix(matrixToPlot = binMatList$Treatment,binSize = opt$bin,axStep = 5000000 ,chromA = "Chr5",startA = 0,endA = 0,chromB = "Chr5",startB = 0,endB = 0)
    dev.off()
}

## Relative difference plot ----

pdf(file = paste(opt$out,"_WGRelDiffMap.pdf",sep=""),height = 19,width = 19)
Genome.f.plot.relative.difference(
    binMatList[[2]],
    binMatList[[1]],
    binSize = opt$bin,
    axStep = 10000000 ,
    chromA = "ALL",
    startA = 0,
    endA = 0,
    chromB = "ALL",
    startB = 0,
    endB = 0,filterThreshold = .99
) -> RelDiffMatrix
dev.off()

if (opt$CHR) {
    pdf(file = paste(opt$out, "_PerChrRelDiffMap.pdf", sep = ""),width = 19,height = 19)
    Chr.f.plot.relative.difference(binMatList[[2]],binMatList[[1]],binSize = opt$bin,axStep = 5000000 ,chromA = "Chr1",startA = 0,endA = 0,chromB = "Chr1",startB = 0,endB = 0)
    Chr.f.plot.relative.difference(binMatList[[2]],binMatList[[1]],binSize = opt$bin,axStep = 5000000 ,chromA = "Chr2",startA = 0,endA = 0,chromB = "Chr2",startB = 0,endB = 0)
    Chr.f.plot.relative.difference(binMatList[[2]],binMatList[[1]],binSize = opt$bin,axStep = 5000000 ,chromA = "Chr3",startA = 0,endA = 0,chromB = "Chr3",startB = 0,endB = 0)
    Chr.f.plot.relative.difference(binMatList[[2]],binMatList[[1]],binSize = opt$bin,axStep = 5000000 ,chromA = "Chr4",startA = 0,endA = 0,chromB = "Chr4",startB = 0,endB = 0)
    Chr.f.plot.relative.difference(binMatList[[2]],binMatList[[1]],binSize = opt$bin,axStep = 5000000 ,chromA = "Chr5",startA = 0,endA = 0,chromB = "Chr5",startB = 0,endB = 0)
    dev.off()
}

## Pre Whole Genome Compare data ----
MA<-binMatList[[1]]
MB<-binMatList[[2]]
diag(MA)<-NA
diag(MB)<-NA
MA[lower.tri(MA)] <- 0
MB[upper.tri(MB)] <- 0
MM<-MA+MB
## Whole Genome Compare ----
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

pdf(file = paste(opt$out,"_WG_CompareMap.pdf",sep=""),width = 19,height = 19)
Genome.f.plot.XY.matrix(
    matrixToPlot = MM,
    binSize = opt$bin,
    axStep = 10000000 ,
    chromA = "ALL",
    startA = 0,
    endA = 0,
    chromB = "ALL",
    startB = 0,
    endB = 0,drawGrid = T
)-> RawIntControl
dev.off()

if (opt$CHR) {
    pdf(file = paste(opt$out, "_PerChrCompareMap.pdf", sep = ""),width = 19,height = 19)
    Chrf.plot.XY.matrix(matrixToPlot = MM,binSize = opt$bin,axStep = 5000000 ,chromA = "Chr1",startA = 0,endA = 0,chromB = "Chr1",startB = 0,endB = 0)
    Chrf.plot.XY.matrix(matrixToPlot = MM,binSize = opt$bin,axStep = 5000000 ,chromA = "Chr2",startA = 0,endA = 0,chromB = "Chr2",startB = 0,endB = 0)
    Chrf.plot.XY.matrix(matrixToPlot = MM,binSize = opt$bin,axStep = 5000000 ,chromA = "Chr3",startA = 0,endA = 0,chromB = "Chr3",startB = 0,endB = 0)
    Chrf.plot.XY.matrix(matrixToPlot = MM,binSize = opt$bin,axStep = 5000000 ,chromA = "Chr4",startA = 0,endA = 0,chromB = "Chr4",startB = 0,endB = 0)
    Chrf.plot.XY.matrix(matrixToPlot = MM,binSize = opt$bin,axStep = 5000000 ,chromA = "Chr5",startA = 0,endA = 0,chromB = "Chr5",startB = 0,endB = 0)
    dev.off()
}

cat(yellow("Finished!"),sep = "\n")
