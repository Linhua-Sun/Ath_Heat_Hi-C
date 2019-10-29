#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(crayon))

option_list = list(
    make_option(
        c("-i", "--input"),
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

if (is.null(opt$input)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

cat(yellow(paste("The input file is ", opt$input,  sep = "")),sep = "\n")
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
    input = opt$input
)

# load all samples into a list, access an individual HiC matrix with binMatList[[sampleName]] or binMatList$sampleName

binMat<-f.load.one.sample(getwd(),files =opt$input,repetitions = 50,binSize = opt$bin )

#binMatList <- f.load.samples(getwd(), sampleList, binSize, 50) ## 5959*5959 MATRIX

write.fst(as.data.frame(binMat),  path = paste(opt$out,"_",opt$input,"_binMat.fst",sep = ""),compress = 100)

## input and Treatment Raw Normalized Heatmap ----

pdf(file = paste(opt$out,"_",opt$input,"_WholeGenomeRawHeatmap.pdf",sep=""),width = 19,height = 19)
Mf.plot.XY.matrix(
    matrixToPlot=binMat,
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
)-> RawIntinput
dev.off()

