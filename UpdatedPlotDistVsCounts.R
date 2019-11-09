#!/usr/bin/env Rscript

## ----
## Linhua Sun
## 20191108｜20191109
## NoScalehicPlotDistVsCounts ==>> ThisScriptRePlot
## ToDo：Use Color random； now the 6 colors are dertermined mannuly.

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
package_list <- c("data.table","ggplot2","scales","ggpmisc","ggsci","ggpubr","sitools","broom")
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

## Import R environments and Define R functions ----

source("/data1/linhua/QIANLAB/R_Functions/theme_linhua.R")

theme_set(theme_linhua(
    base_size = 15,
    legend = "right",
    x.text.angle = 0
))

theme_update(
    #panel.grid.minor.x = element_blank(),
    #panel.grid.minor.y = element_blank(),
    #panel.grid.major.x = element_blank(),
    #panel.grid.major.y = element_blank(),
    #axis.text.x = element_blank(),
    #axis.title.x = element_blank(),
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(0.5, "lines"),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    #aspect.ratio = 1/0.618034,
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0, size = 18),
    legend.title = element_blank(),
    axis.text.y = element_text(hjust = 0),
    axis.text.x = element_text(vjust = 1)
)
base_breaks <- function(n = 10) {
    function(x) {
        axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
    }
}

f2si <-function (number, unit = "")
{
    sifactor <- c(1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-09, 1e-06,
                  0.001, 1, 1000, 1e+06, 1e+09, 1e+12, 1e+15, 1e+18, 1e+21,
                  1e+24)
    pre <- c("y","z","a","f","p","n","u","m","","k","M","G","T","P","E","Z","Y")
    absolutenumber <- number * sign(number)
    ix <- findInterval(absolutenumber, sifactor)
    if (length(ix) > 0) {
        sistring <- paste(number/sifactor[ix], pre[ix], sep = "",
                          unit = unit)
    }
    else {
        sistring <- as.character(number)
    }
    return(sistring)
}
## Import data table ----

opt<-NULL
opt$input<-"ArmsKR_300KToSmKr002K_PerChr.txt"
opt$out<-"PP"
opt$single<-"single"

CountsVsDist<-fread(opt$input,header = F)[V2 !="Matrix"][]

if (unique(CountsVsDist$V3 == "all")==T) {
    opt$single <- "single"
    
} else {
    opt$single <- "NN"
}

CountsVsDist$V5<-as.numeric(CountsVsDist[V2 !="Matrix"]$V5)
CountsVsDist$V4<-as.numeric(CountsVsDist[V2 !="Matrix"]$V4)
colnames(CountsVsDist)<-c("ID", "Matrix", "Chromosome", "Distance", "Contacts")

if (opt$single=="single") {
    
    CountsVsDist[!(Distance %in% c(0))] %>% group_by(Matrix)  %>% do(tidy(lm( log10(Contacts) ~ log10(Distance) , data = .))) %>% as.data.table() %$%  .[term=="log10(Distance)"][,c("Matrix","estimate")] %>% fwrite(paste(opt$out, "_Intercept_AllChr.csv", sep = ""))
    
    ## Single or Per Chromosomes ----
    OutName <- paste(opt$out, "_DecayAllChr.pdf", sep = "")
    #pdf(file = OutName,
    #    height = 6,
    #    width = 8)
    
    formula <- y ~ x
    p2 <-
        ggplot(CountsVsDist[!(Distance %in% c(0))], aes(x = Distance, y = Contacts, color = Matrix)) +
        #geom_smooth(
        #    method = "lm",
        #    formula = formula,
        #    se = F,
        #    size = 0.5
        #) +
        #geom_point(alpha = 10 / 10, size = 2) +
        geom_line(alpha = 1, size = 1) +
        #scale_shape_manual(values = c(0, 1, 2, 5, 6, 21, 22, 23, 24, 25)) +
        stat_poly_eq(
            #aes(label = paste(..rr.label.., ..eq.label.., sep = "~~~~")),
            aes(label = ..eq.label..),
            formula = formula,
            size = 4,
            rr.digits = 3,
            parse = TRUE,
            label.y.npc = "top",
            label.x.npc = "right"
        ) +
        scale_x_log10(name = "Genomic distance (bp)",
                      #breaks = trans_breaks("log10", function(x) 10 ^ x),
                      #breaks=c(10000,100000,1000000,10000000),
                      breaks = base_breaks(n = 4),
                      labels = f2si) +
        scale_y_log10(
            name = "Corrected contact counts",
            breaks = trans_breaks("log10", function(x)
                10 ^ x, n = 5),
            labels = trans_format("log10", math_format(10 ^ .x))
        ) +
        #scale_color_brewer(type = "diverging", palette = "Paired")
        scale_color_manual(
            values =    c(pal_aaas(alpha = 0.7)(9)[5],
                          pal_aaas(alpha = 0.7)(9)[1],
                          pal_aaas(alpha = 0.7)(9)[3],
                          pal_aaas(alpha = 0.7)(9)[6],
                          pal_aaas(alpha = 0.7)(9)[8],
                          pal_aaas(alpha = 0.7)(9)[2]
            )
        )
    p2 + annotation_logticks()
    ggsave(filename = OutName,
           height = 6,
           width = 8)
    
    
} else {
    ## facet version ggplot2 plot ----
    
    CountsVsDist[!(Distance %in% c(0))] %>% group_by(Matrix,Chromosome)  %>% do(tidy(lm( log10(Contacts) ~ log10(Distance) , data = .))) %>% as.data.table() %$%  .[term=="log10(Distance)"][,c("Matrix","Chromosome","estimate")] %>% fwrite(paste(opt$out, "_Intercept_PerChr.csv", sep = ""))
    
    OutName <- paste(opt$out, "_DecayPerChr.pdf", sep = "")
    
    #pdf(file = OutName,
    #    height = 12,
    #    width = 24)
    formula <- y ~ x
    p2 <-
        ggplot(CountsVsDist[!(Distance %in% c(0))], aes(x = Distance, y = Contacts, color = Matrix)) +
        #geom_smooth(
        #    method = "lm",
        #    formula = formula,
        #    se = F,
        #    size = 0.5
        #) +
        #geom_point(alpha = 10 / 10, size = 1) +
        geom_line(alpha = 10 / 10, size = 1) +
        #scale_shape_manual(values = c(0, 1, 2, 5, 6, 21, 22, 23, 24, 25)) +
        stat_poly_eq(
            #aes(label = paste(..rr.label.., ..eq.label.., sep = "~~~~")),
            aes(label = ..eq.label..),
            formula = formula,
            size = 4,
            rr.digits = 3,
            parse = TRUE,
            label.y.npc = "top",
            label.x.npc = "right"
        ) +
        scale_x_log10(name = "Genomic distance (bp)",
                      #breaks = trans_breaks("log10", function(x) 10 ^ x),
                      breaks = base_breaks(5),
                      labels = f2si) +
        scale_y_log10(
            name = "Corrected contact counts",
            breaks = trans_breaks("log10", function(x)
                10 ^ x),
            labels = trans_format("log10", math_format(10 ^ .x))
        ) +
        #scale_color_ucscgb(alpha = 1)
        scale_color_manual(
            values =    c(pal_aaas(alpha = 0.7)(9)[5],
                          pal_aaas(alpha = 0.7)(9)[1],
                          pal_aaas(alpha = 0.7)(9)[3],
                          pal_aaas(alpha = 0.7)(9)[6],
                          pal_aaas(alpha = 0.7)(9)[8],
                          pal_aaas(alpha = 0.7)(9)[2]
            )
        )
    p2 + annotation_logticks() -> PP
    facet(p = PP, facet.by = "Chromosome")
    #graphics.off()
    ggsave(filename = OutName,height = 12,width = 24)
}