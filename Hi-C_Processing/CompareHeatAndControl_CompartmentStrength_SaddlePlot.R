## Aim: Chromatin Compartment Strength Analaysis
## Linhua Sun
## 

library(ggpubr)
library(ggsci)
library(GENOVA)

## Define new function to return values ----

New_visualise.saddle <-
    function (SBoutList,
              addText = T,
              zlim = c(0.5, 2),
              EVlim = c(-1.5,1.5),
              square = T,
              crossLines = NULL)
        
    {
        library(fields)
        if (is.null(EVlim)) {
            EVlim = zlim
        }
        df_mat = NULL
        df_ev = NULL
        NSAMPLE = length(SBoutList)
        SAMPLENAMES = c()
        for (i in 1:NSAMPLE) {
            SAMPLENAMES = c(SAMPLENAMES, as.character(unique(SBoutList[[i]]$MAT$sample)))
            df_mat = rbind(df_mat, SBoutList[[i]]$MAT)
            df_ev = rbind(df_ev, SBoutList[[i]]$EV)
        }
        if (is.null(crossLines)) {
            if (max(df_mat$Var1) < 10) {
                warning("Number of bins is lower that 10: not showing crossLines")
                crossLines = F
            }
            else {
                crossLines = T
            }
        }
        lay = matrix(0, nrow = 4, ncol = NSAMPLE)
        lay[1, ] = 1:NSAMPLE
        for (i in 1:NSAMPLE) {
            lay[2:4, i] = i + NSAMPLE
        }
        layout(lay)
        par(mar = rep(1, 4), xaxs = "i", yaxs = "i")
        EV0 = c()
        newMat <- NULL
        df_ev$ID <- factor(apply(df_ev[, 3:5], 1, paste0, collapse = "_"))
        for (i in unique(df_ev$ID)) {
            tt <- df_ev[df_ev$ID == i, ]
            ttscaled <- as.data.frame(approx(tt$ev, n = 100))
            ttscaled <- suppressWarnings(cbind(ttscaled, unique(df_ev[df_ev$ID == 
                                                                          i, c("chrom", "arm", "sample", "ID")])))
            ttscaled$x <- 1:100
            newMat <- rbind(newMat, ttscaled)
        }
        toPlotEV <- NULL
        for (s in unique(newMat$sample)) {
            tmp <- newMat[newMat$sample == s, ]
            loesje <- as.data.frame(loess.smooth(x = tmp$x, y = tmp$y))
            loesje$sample = s
            toPlotEV <- rbind(toPlotEV, loesje)
        }
        toPlotEV <- toPlotEV[, c(3, 1, 2)]
        for (S in SAMPLENAMES) {
            tmp = setNames(toPlotEV[toPlotEV[, 1] == S, 2:3], c("x", 
                                                                "y"))
            tmp.bk = tmp
            crossPoint = approx(y = tmp.bk$x, x = tmp.bk$y, xout = 0)$y
            EV0 = c(EV0, crossPoint/max(toPlotEV[, 2]))
            X <- tmp.bk$x
            y.low <- rep(0, length(tmp.bk$x))
            y.high <- tmp.bk$y
            plot(X, y.high, type = "n", ylim = EVlim, axes = F, main = S, 
                 ylab = "", xlab = "")
            lines(X, y.low, col = "black")
            lines(X, y.high, col = "black")
            polygon(c(X, rev(X)), c(y.high, rev(y.low)), col = "black", 
                    border = NA)
            box()
            if (crossLines) {
                abline(v = crossPoint, lty = 2)
            }
        }
        if (square) {
            par(pty = "s")
        }
        PAL = colorRampPalette(rev(c("#B2182B", "white", "#2166AC")))
        for (Si in 1:NSAMPLE) {
            S = SAMPLENAMES[Si]
            tmp = df_mat[df_mat$sample == S, ]
            tmp = setNames(aggregate(tmp$value, by = list(tmp$Var1, 
                                                          tmp$Var2), mean, na.rm = T), c("x", "y", "z"))
            M = matrix(0, nrow = max(tmp$x), ncol = max(tmp$y))
            M[cbind(tmp$x, tmp$y)] = tmp$z
            M[M < zlim[1]] = zlim[1]
            M[M > zlim[2]] = zlim[2]
            image(M, zlim = zlim, col = PAL(50), axes = F, breaks = 2^(seq(log2(zlim[1]), 
                                                                           log2(zlim[2]), length.out = 51)))
            box()
            if (addText) {
                text(x = 0.05, y = 0.95, label = "AB")
                text(y = 0.05, x = 0.95, label = "BA")
                text(y = 0.05, x = 0.05, label = "BB")
                text(y = 0.95, x = 0.95, label = "AA")
            }
            if (crossLines) {
                abline(v = (EV0[Si]) - (0.5/max(tmp$y)), lty = 2)
                abline(h = ((EV0[Si]) - (0.5/max(tmp$y))), lty = 2)
            }
        }
        if (square) {
            par(pty = "m")
        }
        return(df_mat)
    }

Edited.visualise.compartmentStrength <- function (SBoutList, showInteractions = F) 
{
    require(ggplot2)
    strengthDF = data.frame()
    interactionDF = data.frame()
    namesVector <- c()
    for (i in 1:length(SBoutList)) {
        dat = SBoutList[[i]]
        namesVector <- c(namesVector, unique(dat$MAT$sample))
        dat$MAT$CC = "XX"
        MAXbin = max(dat$MAT$Var1)
        binsTOse = floor(MAXbin * 0.2)
        binsTOse = max(1, binsTOse)
        dat$MAT$unLog = 2^dat$MAT$value
        dat$MAT[dat$MAT$Var1 <= binsTOse & dat$MAT$Var2 <= binsTOse, 
                "CC"] = "BB"
        dat$MAT[dat$MAT$Var2 <= binsTOse & dat$MAT$Var1 >= MAXbin - 
                    binsTOse + 1, "CC"] = "AB"
        dat$MAT[dat$MAT$Var1 >= MAXbin - binsTOse + 1 & dat$MAT$Var2 >= 
                    MAXbin - binsTOse + 1, "CC"] = "AA"
        dat$MAT = dat$MAT[dat$MAT$CC != "XX", ]
        tmp = dplyr::summarise(dplyr::group_by(dat$MAT, color, 
                                               sample, chrom, arm, CC), score = mean(unLog))
        interactionDF = rbind(interactionDF, as.data.frame(tmp))
        for (S in unique(tmp$sample)) {
            meta = unique(tmp[tmp$sample == S, c("chrom", "arm")])
            for (C in unique(meta$chrom)) {
                for (A in unique(unname(unlist(meta[meta$chrom == 
                                                    C, "arm"])))) {
                    tmpi = tmp[tmp$sample == S & tmp$chrom == 
                                   C & tmp$arm == A, ]
                    strength = log(tmpi[tmpi$CC == "AA", "score"] * 
                                       tmpi[tmpi$CC == "BB", "score"]/tmpi[tmpi$CC == 
                                                                               "AB", "score"]^2)
                    strengthDF = rbind(strengthDF, data.frame(S, 
                                                              C, A, strength, unique(tmp[tmp$sample == 
                                                                                             S, "color"])))
                }
            }
        }
    }
    RANGE = range(strengthDF$score, na.rm = T)
    RANGE[1] = floor(RANGE[1]) * 0.95
    RANGE[2] = ceiling(RANGE[2]) * 1.05
    if (showInteractions != TRUE) {
        coltmp <- col2rgb(levels(as.factor(strengthDF$color)), 
                          alpha = T)/255
        coltmp[4, ] <- coltmp[4, ] * 0.85
        cols <- rgb(red = coltmp[1, ], green = coltmp[2, ], 
                    blue = coltmp[3, ], alpha = coltmp[4, ])
        boxplot(split(strengthDF$score, strengthDF$S), ylim = RANGE, 
                col = cols, ylab = "compartment strength")
    }
    else {
        interactionDF$sample <- factor(interactionDF$sample, 
                                       levels = namesVector)
        P <- ggplot2::ggplot(interactionDF, ggplot2::aes(x = sample, 
                                                         y = score, fill = color)) + ggplot2::facet_wrap("CC") + 
            ggplot2::geom_hline(yintercept = 1, lty = 3) + ggplot2::geom_boxplot() + 
            ggplot2::labs(y = "Average O/E") + GENOVA_THEME() + 
            ggplot2::scale_fill_identity() + ggplot2::guides(fill = F) + 
            ggplot2::theme(axis.title.x = ggplot2::element_blank())
        print(P)
    }
    return(list(strengthDF, interactionDF))
}
## Import data ----

setwd("/data1/linhua/QIANLAB/PROJECT/hic/data")
BED20K="/data1/linhua/QIANLAB/PROJECT/hic/Matrix_Heat/Col22_total_sample/raw/20000/Col22_total_sample_20000_abs.bed"
Ath_centromeres = read.delim('/data1/linhua/QIANLAB/PROJECT/hic/Ath_Centromere.bed', sep = '\t', h = F, stringsAsFactors = F)

H3K4me3<-"/data1/linhua/QIANLAB/PROJECT/hic/Hi-C-2019/HM_H3K4me3_Ath_rosette_leaf_normal.bed"
H3K4me3Peaks = read.delim(H3K4me3, h = F)

Control_20Kb <-
    construct.experiment(
        ignore.checks = T,
        signalPath = '/data1/linhua/QIANLAB/PROJECT/hic/HeatProjectHi-C-Results-ToBeSubmitted/Control_20000_iced.matrix',
        indicesPath = BED20K,
        centromeres = Ath_centromeres,
        name = "Control",
        color = "black"
    )

Heat_20Kb <- construct.experiment(
    ignore.checks = T,
    signalPath = '/data1/linhua/QIANLAB/PROJECT/hic/HeatProjectHi-C-Results-ToBeSubmitted/Heat_20000_iced.matrix',
    indicesPath = BED20K,
    centromeres = Ath_centromeres,
    name = "Heat",
    color = "red"
)



saddleControl_20Kb = saddle(exp = Control_20Kb,chip = H3K4me3Peaks,chromsToUse = paste0('Chr', 1:5),nBins = 50)
saddleHeat_20Kb = saddle(exp = Heat_20Kb,chip = H3K4me3Peaks,chromsToUse = paste0('Chr', 1:5),nBins = 50)

## saddle plot ----

pdf(file = "PooledHeatVSControlsaddle.pdf")
df_mat <- New_visualise.saddle(SBoutList = list(saddleControl_20Kb,saddleHeat_20Kb),
                               crossLines = T,
                               addText = T)
dev.off()

## Plot By Chromosome arms ----
library(ggsci)

PP<-Edited.visualise.compartmentStrength(list(saddleControl_20Kb,saddleHeat_20Kb))

LL<-PP[[1]]
II<-PP[[2]]

LL$GG<-paste(LL$C,LL$A,sep = "")
II$GG<-paste(II$chrom,II$arm,sep = "")

pdf(file = "HeatVSControlHomo.pdf",height = 6,width = 4)
LL %>%
    ggplot(aes(x=S, y=score,fill=S)) +
    geom_boxplot(width=0.5)+
    geom_line(aes(x=S, y=score,group=GG),color="grey")+
    geom_point(size=2,color="black",fill="grey",shape=21,alpha=1)+
    scale_fill_manual(values = pal_aaas()(9)[3:2])+
    ylab("Compartment strength")+xlab("Conditions")+
    #geom_hline(yintercept = 1,color="black",size=0.3,linetype="dashed")+
    #geom_hline(yintercept = 0,color="black",size=0.3,linetype="dashed")+
    theme_classic(base_size = 15)+theme(aspect.ratio = 2,legend.position="none")+border()+border()
dev.off()

pdf(file = "HeatVSControlHetero.pdf",height = 6,width = 12)
II %>%
    ggplot(aes(x=sample, y=score,fill=sample)) +
    geom_boxplot(width=0.5)+
    geom_line(aes(x=sample, y=score,group=GG),color="grey")+
    geom_point(size=2,color="black",fill="grey",shape=21,alpha=1)+
    scale_fill_manual(values = pal_aaas()(9)[3:2])+
    ylab("Average O/E")+xlab("Conditions")+
    #geom_hline(yintercept = 1,color="black",size=0.3,linetype="dashed")+
    #geom_hline(yintercept = 0,color="black",size=0.3,linetype="dashed")+
    theme_classic(base_size = 15)+theme(aspect.ratio = 2,legend.position="none")+border() -> PLOT
PLOT  %>%  facet(facet.by = "CC")
dev.off()
