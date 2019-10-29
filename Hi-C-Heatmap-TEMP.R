#!/usr/bin/env Rscript
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

source("/data1/linhua/QIANLAB/PROJECT/hic/FPC/BIN100K_Col_Matrix/BU/Modified_Plot_Heatmap_Base.R")

## Input and Output ----
pathToScripts <- "/data1/linhua/QIANLAB/PROJECT/hic/MWSchmid_test/Rscripts/"
f.source.organism.specific.code(file.path(pathToScripts, "HiCdat-A-thaliana-TAIR10.R"))
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

pdf(file = paste(opt$out,"_",opt$Control,"_WholeGenomeRawHeatmap.pdf",sep=""),width = 19,height = 19)
Mf.plot.XY.matrix(
    matrixToPlot=binMatList$Control,
    binSize =opt$bin,
    axStep =opt$bin ,
    chromA = "ALL",
    startA = 0,
    endA = 0,
    chromB = "ALL",
    startB = 0,
    endB = 0,
    useLog = TRUE,
    drawGrid = TRUE,
    doNorm = FALSE,
    doCor = F,
    useSplineInterPol = TRUE
)-> RawIntControl
dev.off()

pdf(file = paste(opt$out,"_",opt$Treatment,"_WholeGenomeRawHeatmap.pdf",sep=""),width = 19,height = 19)
Mf.plot.XY.matrix(
    matrixToPlot=binMatList$Treatment,
    binSize =opt$bin,
    axStep =opt$bin ,
    chromA = "ALL",
    startA = 0,
    endA = 0,
    chromB = "ALL",
    startB = 0,
    endB = 0,
    useLog = TRUE,
    drawGrid = TRUE,
    doNorm = FALSE,
    doCor = F,
    useSplineInterPol = TRUE
) -> RawIntTreatment
dev.off()

## Relative difference plot ----
pdf(file = paste(opt$out,"_WholeGenomeRelDifferenceHeatmap.pdf",sep=""),height = 19,width = 19)
Mf.plot.relative.difference(
    binMatList[[2]],
    binMatList[[1]],
    binSize,
    pathToTutorial,
    "mQIAN_Diff_relative_Col_vs_Cold_100K"
)-> RelDiffMatrix
dev.off()

## Pre Whole Genome Compare data ----
MA<-binMatList[[1]]
MB<-binMatList[[2]]
diag(MA)<-NA
diag(MB)<-NA
MA[lower.tri(MA)] <- 0
MB[upper.tri(MB)] <- 0
MM<-MA+MB
## Whole Genome Compare ----
pdf(file = paste(opt$out,"_WholeGenomeCompareHeatmap.pdf",sep=""),width = 19,height = 19)
Mf.plot.XY.matrix(
    matrixToPlot=MM,
    binSize =opt$bin,
    axStep =opt$bin ,
    chromA = "ALL",
    startA = 0,
    endA = 0,
    chromB = "ALL",
    startB = 0,
    endB = 0,
    useLog = TRUE,
    drawGrid = TRUE,
    doNorm = FALSE,
    doCor = F,
    useSplineInterPol = TRUE
)-> WholeGenomeCompareMatrix
dev.off()

## ----
cat(yellow("Finished!"),sep = "\n")
