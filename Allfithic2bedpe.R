#!/usr/bin/env Rscript
##############################################################################################################

## Aim in convert fithic results into bedpe for IGV
## Linhua Sun 
## 20191008

library(optparse)
suppressWarnings(suppressMessages(library("crayon",character.only = TRUE,quietly = TRUE,warn.conflicts = FALSE)))
option_list = list(
    make_option(
        c("-i", "--input"),
        type = "character",
        default = NULL,
        help = "FitHiC results file name",
        metavar = "character"
    ),
    make_option(
        c("-o", "--out"),
        type = "character",
        default = "outputprefix",
        help = "outputprefix [default= %default]",
        metavar = "character"
    )
)

opt_parser = OptionParser(option_list = option_list)

opt = parse_args(opt_parser)

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop(red("At least 3 argument must be supplied: Rscriipt ......"), call.=FALSE)
}

cat(yellow(paste("The input file is ", opt$input,  sep = "")),sep = "\n")
cat(yellow(paste("The output prefix is ", opt$out, sep = "")),sep = "\n")

###################################################################

# See whether these packages exist on comp. If not, install.
package_list <- c("data.table")

for(p in package_list) {
    if (!suppressWarnings(suppressMessages(require(
        p,
        character.only = TRUE,
        quietly = TRUE,
        warn.conflicts = FALSE
    )))) {
        install.packages(p, repos = "http://cran.r-project.org")
        suppressWarnings(suppressMessages(library(
            p,
            character.only = TRUE,
            quietly = TRUE,
            warn.conflicts = FALSE
        )))
    }
}

## use pacman as alternative to import a series R packages
###################################################################

ExtractLoopsFitHiC <- function(Resfile,OutName) {
    TT <-fread(Resfile)
    
    #PP <- TT[, c("chr1", "fragmentMid1", "fragmentMid2","q-value")]
    #PP$name <- paste("stemloop", seq(1, nrow(PP)), sep = "")
    #write("track graphType=arc",file =  paste(OutName,"_arc_01.bed",sep=""))
    #fwrite(PP[`q-value` < 0.01][,c("chr1", "fragmentMid1", "fragmentMid2","name")], paste(OutName,"_arc_01.bed",sep=""), sep = "\t", col.names =  F,append = T)
    #system(paste("bgzip -f ",paste(OutName,"_arc_01.bed",sep="")))
    ##system(paste("tabix -f -p bed ",paste(OutName,"_arc_01.bed.gz",sep="")))
    
    TT$start1<-TT$fragmentMid1-1
    TT$end1<-TT$fragmentMid1+1
    TT$start2<-TT$fragmentMid2-1
    TT$end2<-TT$fragmentMid2+1
    TT$name<-"."
    #TT[TT$contactCount>1000]$contactCount<-1000
    TT$score<-(-log10(TT$`q-value`))
    
    TT[`q-value`<0.01][,c("chr1","start1","end1","chr2","start2","end2","score")] %>% fwrite(file = paste(OutName,"_Loop_01.bedpe",sep=""),sep = "\t",col.names = F)
    system(paste("bgzip -f ",paste(OutName,"_Loop_01.bedpe",sep="") ))
    system(paste("tabix -f -p bed ",paste(OutName,"_Loop_01.bedpe.gz",sep="")))
  
    TT[,c("chr1","start1","end1","chr2","start2","end2","score")] %>% fwrite(file = paste(OutName,"_Loop_All.bedpe",sep=""),sep = "\t",col.names = F)
    system(paste("bgzip -f ",paste(OutName,"_Loop_All.bedpe",sep="") ))
    system(paste("tabix -f -p bed ",paste(OutName,"_Loop_All.bedpe.gz",sep="")))
    
}

ExtractLoopsFitHiC(opt$input,OutName = opt$out)