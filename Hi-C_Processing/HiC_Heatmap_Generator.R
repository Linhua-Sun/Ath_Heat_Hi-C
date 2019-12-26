library(fields)

Genome.f.plot.relative.difference <-
    function (dataMatrixA,
              dataMatrixB,
              binSize,
              axStep,
              chromA = "ALL",
              startA = 0,
              endA = 0,
              chromB = "ALL",
              startB = 0,
              endB = 0,
              filterZero = TRUE,
              filterThreshold = 0.95,
              drawGrid = T,
              randomizeDiff = FALSE)

    {
        seList <- f.get.se.list(binSize)
        chromSizes <- f.get.chrom.sizes()
        if (randomizeDiff) {
            rawDiff <- f.internal.randomize.matrix(dataMatrixA -
                                                       dataMatrixB, binSize)
        } else {
            rawDiff <- dataMatrixA - dataMatrixB
        }
        
        if (filterZero) {
            toKeep <- intersect(
                f.internal.find.non.zero.indices(dataMatrixA,
                                                 filterThreshold),
                f.internal.find.non.zero.indices(dataMatrixB,
                                                 filterThreshold)
            )
        } else {
            toKeep <- 1:nrow(rawDiff)
        }
        
        if (filterZero) {
            toRemove <- setdiff(1:nrow(rawDiff), toKeep)
            rawDiff[toRemove,] <- 0
            rawDiff[, toRemove] <- 0
        }
        ## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if (startA < axStep) {
            startA <- 0
        }
        
        if (startB < axStep) {
            startB <- 0
        }
        
        if (chromA == "ALL") {
            sa <- 1
            ea <- nrow(dataMatrixA)
            xAxAt <-
                as.vector(unlist(
                    f.internal.axis.maker.on.index(binSize,
                                                   axStep, seList)
                ))
            xAxLab <-
                as.vector(unlist(f.internal.axis.maker(binSize,
                                                       axStep, seList))) /
                1e+06
            xGrid <- do.call("rbind", seList)[, 1]
        } else {
            if (endA > chromSizes[[chromA]]) {
                endA <- chromSizes[[chromA]] ## Quality control of the value of endA
            }
            if (endA == 0) {
                endA <- chromSizes[[chromA]]
            }
            sa <- f.translate.chrom.pos.to.index(chromA, startA,
                                                 seList, binSize)
            ea <-
                f.translate.chrom.pos.to.index(chromA, endA, seList,
                                               binSize)
            xAxAt <- seq(0, (ea - sa), by = axStep / binSize)
            xAxLab <- seq(startA, endA, by = axStep) / 1e+06
            xGrid <- xAxAt
        }
        if (chromB == "ALL") {
            sb <- 1
            eb <- ncol(dataMatrixA)
            yAxAt <-
                as.vector(unlist(
                    f.internal.axis.maker.on.index(binSize,
                                                   axStep, seList)
                ))
            yAxLab <-
                as.vector(unlist(f.internal.axis.maker(binSize,
                                                       axStep, seList))) /
                1e+06
            yGrid <- do.call("rbind", seList)[, 1]
        } else {
            if (endB > chromSizes[[chromB]]) {
                endB <- chromSizes[[chromB]]
            }
            if (endB == 0) {
                endB <- chromSizes[[chromB]]
            }
            sb <- f.translate.chrom.pos.to.index(chromB, startB,
                                                 seList, binSize)
            eb <-
                f.translate.chrom.pos.to.index(chromB, endB, seList,
                                               binSize)
            yAxAt <- seq(0, (eb - sb), by = axStep / binSize)
            yAxLab <- seq(startB, endB, by = axStep) / 1e+06
            yGrid <- yAxAt
        }
        
        ## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        
        print(paste(sa, ea, sb, eb, sep = " "))
        
        normDiff <- rawDiff / ((dataMatrixA + dataMatrixB) / 2)
        normDiff[is.na(normDiff)] <- 0
        
        #normDiff[normDiff>1]<-1
        #normDiff[normDiff<(-1)]<- (-1)
        
        normDiff <-  normDiff[sa:ea, sb:eb]
        
        par(oma = c(5, 5, 2, 2), mar = c(5, 5, 2,2))
        image.plot(
            1:nrow(normDiff),
            1:ncol(normDiff),
            zlim = c(-1, 1),
            normDiff,
            col = colorRampPalette(c(
                "green", "blue", "white", "red", "orange"
            ))(64),
            useRaster = TRUE,
            axis.args = list(cex.axis = 3),
            yaxt = "n",
            xaxt = "n",
            xlab = "",
            ylab = ""
        )
        
        xDT <- data.table(xAxAt, xAxLab)
        yDT <- data.table(yAxAt, yAxLab)
        
        axis(1,at = xDT[xAxLab != 30]$xAxAt,labels = xDT[xAxLab != 30]$xAxLab,outer = FALSE,line = 0,lwd = 2,cex.axis = 2,las = 1)
        #axis(2, at = axAt, labels = axLab, outer = FALSE, line = 2,lwd = 2, cex.axis = 2, las = 1)
        
        if (drawGrid) {
            abline(
                h = yGrid,
                v = xGrid,
                lwd = 2,
                col = "black"
            )
        }
        
        return(normDiff)
        
        
        
    }

Chr.f.plot.relative.difference <-
    function (dataMatrixA,
              dataMatrixB,
              binSize,
              axStep,
              chromA = "ALL",
              startA = 0,
              endA = 0,
              chromB = "ALL",
              startB = 0,
              endB = 0,
              filterZero = TRUE,
              filterThreshold = 0.95,
              drawGrid = T,
              randomizeDiff = FALSE)

    {
        seList <- f.get.se.list(binSize)
        chromSizes <- f.get.chrom.sizes()
        if (randomizeDiff) {
            rawDiff <- f.internal.randomize.matrix(dataMatrixA -
                                                       dataMatrixB, binSize)
        } else {
            rawDiff <- dataMatrixA - dataMatrixB
        }
        
        if (filterZero) {
            toKeep <- intersect(
                f.internal.find.non.zero.indices(dataMatrixA,
                                                 filterThreshold),
                f.internal.find.non.zero.indices(dataMatrixB,
                                                 filterThreshold)
            )
        } else {
            toKeep <- 1:nrow(rawDiff)
        }
        
        if (filterZero) {
            toRemove <- setdiff(1:nrow(rawDiff), toKeep)
            rawDiff[toRemove,] <- 0
            rawDiff[, toRemove] <- 0
        }
        ## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if (startA < axStep) {
            startA <- 0
        }
        
        if (startB < axStep) {
            startB <- 0
        }
        
        if (chromA == "ALL") {
            sa <- 1
            ea <- nrow(dataMatrixA)
            xAxAt <-
                as.vector(unlist(
                    f.internal.axis.maker.on.index(binSize,
                                                   axStep, seList)
                ))
            xAxLab <-
                as.vector(unlist(f.internal.axis.maker(binSize,
                                                       axStep, seList))) /
                1e+06
            xGrid <- do.call("rbind", seList)[, 1]
        } else {
            if (endA > chromSizes[[chromA]]) {
                endA <- chromSizes[[chromA]] ## Quality control of the value of endA
            }
            if (endA == 0) {
                endA <- chromSizes[[chromA]]
            }
            sa <- f.translate.chrom.pos.to.index(chromA, startA,
                                                 seList, binSize)
            ea <-
                f.translate.chrom.pos.to.index(chromA, endA, seList,
                                               binSize)
            xAxAt <- seq(0, (ea - sa), by = axStep / binSize)
            xAxLab <- seq(startA, endA, by = axStep) / 1e+06
            xGrid <- xAxAt
        }
        if (chromB == "ALL") {
            sb <- 1
            eb <- ncol(dataMatrixA)
            yAxAt <-
                as.vector(unlist(
                    f.internal.axis.maker.on.index(binSize,
                                                   axStep, seList)
                ))
            yAxLab <-
                as.vector(unlist(f.internal.axis.maker(binSize,
                                                       axStep, seList))) /
                1e+06
            yGrid <- do.call("rbind", seList)[, 1]
        } else {
            if (endB > chromSizes[[chromB]]) {
                endB <- chromSizes[[chromB]]
            }
            if (endB == 0) {
                endB <- chromSizes[[chromB]]
            }
            sb <- f.translate.chrom.pos.to.index(chromB, startB,
                                                 seList, binSize)
            eb <-
                f.translate.chrom.pos.to.index(chromB, endB, seList,
                                               binSize)
            yAxAt <- seq(0, (eb - sb), by = axStep / binSize)
            yAxLab <- seq(startB, endB, by = axStep) / 1e+06
            yGrid <- yAxAt
        }
        
        ## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        
        print(paste(sa, ea, sb, eb, sep = " "))
        
        normDiff <- rawDiff / ((dataMatrixA + dataMatrixB) / 2)
        normDiff[is.na(normDiff)] <- 0
        
        #normDiff[normDiff>1]<-1
        #normDiff[normDiff<(-1)]<- (-1)
        
        normDiff <-  normDiff[sa:ea, sb:eb]
        
        #par(oma = c(5, 5, 3, 3), mar = c(5, 5, 3, 3))
        par(oma = c(5, 5, 2, 2), mar = c(5, 5, 2,2))
        image.plot(
            1:nrow(normDiff),
            1:ncol(normDiff),
            zlim = c(-1, 1),
            normDiff,
            col = colorRampPalette(c(
                "green", "blue", "white", "red", "orange"
            ))(64),
            useRaster = TRUE,
            axis.args = list(cex.axis = 3),
            yaxt = "n",
            xaxt = "n",
            xlab = "",
            ylab = ""
        )
        xDT <- data.table(xAxAt, xAxLab)
        yDT <- data.table(yAxAt, yAxLab)
        #axis(1,at = xDT[xAxLab != 30]$xAxAt,labels = xDT[xAxLab != 30]$xAxLab,outer = FALSE,line = 0,lwd = 2,cex.axis = 2,las = 1)
        axis(1, at = xDT[xAxLab!=30]$xAxAt, labels = paste(xDT[xAxLab!=30]$xAxLab,"Mb",sep = ""), outer = FALSE, line = 0,lwd = 2, cex.axis = 2, las = 1)
        #axis(2, at = axAt, labels = axLab, outer = FALSE, line = 2,lwd = 2, cex.axis = 2, las = 1)
        #if (drawGrid) {
        #    abline(
        #        h = yGrid,
        #        v = xGrid,
        #        lwd = 2,
        #        col = "black"
        #    )
        #}
        return(normDiff)
    }


## ----
Genome.f.plot.XY.matrix <-
    function (matrixToPlot,
              binSize,
              axStep,
              chromA = "ALL",
              startA = 0,
              endA = 0,
              chromB = "ALL",
              startB = 0,
              endB = 0,
              useLog = TRUE,
              drawGrid = FALSE,
              useSplineInterPol = TRUE
) {
    source("/data1/linhua/QIANLAB/PROJECT/hic/MWSchmid_test/Rscripts//HiCdat-A-thaliana-TAIR10.R")
    seList <- f.get.se.list(binSize)
    chromSizes <- f.get.chrom.sizes()
    if (startA < axStep) {
        startA <- 0
    }
    if (startB < axStep) {
        startB <- 0
    }
    if (chromA == "ALL") {
        sa <- 1
        ea <- nrow(matrixToPlot)
        xAxAt <- as.vector(unlist(f.internal.axis.maker.on.index(binSize,
                                                                 axStep, seList)))
        xAxLab <- as.vector(unlist(f.internal.axis.maker(binSize,
                                                         axStep, seList)))/1e+06
        xGrid <- do.call("rbind", seList)[, 1]
    } else {
        if (endA > chromSizes[[chromA]]) {
            endA <- chromSizes[[chromA]] ## Quality control of the value of endA
        } 
        if (endA ==0 ) {
            endA <- chromSizes[[chromA]] 
        } 
        sa <- f.translate.chrom.pos.to.index(chromA, startA,
                                             seList, binSize)
        ea <- f.translate.chrom.pos.to.index(chromA, endA, seList,
                                             binSize)
        xAxAt <- seq(0, (ea - sa), by = axStep/binSize)
        xAxLab <- seq(startA, endA, by = axStep)/1e+06
        xGrid <- xAxAt
    }
    if (chromB == "ALL") {
        sb <- 1
        eb <- ncol(matrixToPlot)
        yAxAt <- as.vector(unlist(f.internal.axis.maker.on.index(binSize,
                                                                 axStep, seList)))
        yAxLab <- as.vector(unlist(f.internal.axis.maker(binSize,
                                                         axStep, seList)))/1e+06
        yGrid <- do.call("rbind", seList)[, 1]
    } else {
        if (endB > chromSizes[[chromB]]) {
            endB <- chromSizes[[chromB]]
        }
        if (endB ==0 ) {
            endB <- chromSizes[[chromB]] 
        } 
        sb <- f.translate.chrom.pos.to.index(chromB, startB,
                                             seList, binSize)
        eb <- f.translate.chrom.pos.to.index(chromB, endB, seList,
                                             binSize)
        yAxAt <- seq(0, (eb - sb), by = axStep/binSize)
        yAxLab <- seq(startB, endB, by = axStep)/1e+06
        yGrid <- yAxAt
    }
    
    matrixToPlot <- log2(matrixToPlot + 1)
    
    print(paste(sa,ea, sb,eb,sep = " "))
    
    matrixToPlot <- matrixToPlot[sa:ea, sb:eb]
    
    par(oma = c(5, 5, 2, 2), mar = c(5, 5, 2,2))
    
    matrixToPlot[matrixToPlot>8]<-8
    diag(matrixToPlot)<-NA
    image.plot(
        1:nrow(matrixToPlot),
        1:ncol(matrixToPlot),
        matrixToPlot,
        zlim = c(0, 8),
        useRaster = T,
        col = colorRampPalette(c(
            "black", "yellow", "orange", "red", "darkred"
        ))(64),
        axis.args = list(cex.axis = 3),
        yaxt = "n",
        xaxt = "n",
        xlab = "",
        ylab = ""
    )
    
    xDT<-data.table(xAxAt,xAxLab)
    yDT<-data.table(yAxAt,yAxLab)
    
    axis(1, at = xDT[xAxLab!=30]$xAxAt, labels = xDT[xAxLab!=30]$xAxLab, outer = FALSE, line = NA,lwd.ticks = 2, cex.axis = 2, las = 1)
    #axis(2, at = yDT[yAxLab!=30]$yAxAt, labels = paste(yDT[yAxLab!=30]$yAxLab,"Mb",sep = ""), outer = FALSE, line = NA,lwd.ticks = 2, cex.axis = 2, las = 1)
    
    if (drawGrid) {abline(h = yGrid, v = xGrid, lwd = 2, col = "white")}
    return(matrixToPlot)
    }

## ----
Chrf.plot.XY.matrix <-
    function (matrixToPlot,
              binSize,
              axStep,
              chromA = "ALL",
              startA = 0,
              endA = 0,
              chromB = "ALL",
              startB = 0,
              endB = 0,
              useLog = TRUE,
              drawGrid = FALSE
    ) {
        source("/data1/linhua/QIANLAB/PROJECT/hic/MWSchmid_test/Rscripts//HiCdat-A-thaliana-TAIR10.R")
        seList <- f.get.se.list(binSize)
        chromSizes <- f.get.chrom.sizes()
        if (startA < axStep) {
            startA <- 0
        }
        if (startB < axStep) {
            startB <- 0
        }
        if (chromA == "ALL") {
            sa <- 1
            ea <- nrow(matrixToPlot)
            xAxAt <- as.vector(unlist(f.internal.axis.maker.on.index(binSize,
                                                                     axStep, seList)))
            xAxLab <- as.vector(unlist(f.internal.axis.maker(binSize,
                                                             axStep, seList)))/1e+06
            xGrid <- do.call("rbind", seList)[, 1]
        } else {
            if (endA > chromSizes[[chromA]]) {
                endA <- chromSizes[[chromA]] ## Quality control of the value of endA
            } 
            if (endA ==0 ) {
                endA <- chromSizes[[chromA]] 
            } 
            sa <- f.translate.chrom.pos.to.index(chromA, startA,
                                                 seList, binSize)
            ea <- f.translate.chrom.pos.to.index(chromA, endA, seList,
                                                 binSize)
            xAxAt <- seq(0, (ea - sa), by = axStep/binSize)
            xAxLab <- seq(startA, endA, by = axStep)/1e+06
            xGrid <- xAxAt
        }
        if (chromB == "ALL") {
            sb <- 1
            eb <- ncol(matrixToPlot)
            yAxAt <- as.vector(unlist(f.internal.axis.maker.on.index(binSize,
                                                                     axStep, seList)))
            yAxLab <- as.vector(unlist(f.internal.axis.maker(binSize,
                                                             axStep, seList)))/1e+06
            yGrid <- do.call("rbind", seList)[, 1]
        } else {
            if (endB > chromSizes[[chromB]]) {
                endB <- chromSizes[[chromB]]
            }
            if (endB ==0 ) {
                endB <- chromSizes[[chromB]] 
            } 
            sb <- f.translate.chrom.pos.to.index(chromB, startB,
                                                 seList, binSize)
            eb <- f.translate.chrom.pos.to.index(chromB, endB, seList,
                                                 binSize)
            yAxAt <- seq(0, (eb - sb), by = axStep/binSize)
            yAxLab <- seq(startB, endB, by = axStep)/1e+06
            yGrid <- yAxAt
        }
        
        matrixToPlot <- log2(matrixToPlot + 1)
        
        print(paste(sa,ea, sb,eb,sep = " "))
        
        matrixToPlot <- matrixToPlot[sa:ea, sb:eb]
        
        par(oma = c(5, 5, 2, 2), mar = c(5, 5, 2,2))
        
        matrixToPlot[matrixToPlot>8]<-8
        diag(matrixToPlot)<-NA
        image.plot(
            1:nrow(matrixToPlot),
            1:ncol(matrixToPlot),
            matrixToPlot,
            zlim = c(0, 8),
            useRaster = T,
            col = colorRampPalette(c(
                "black", "yellow", "orange", "red", "darkred"
            ))(64),
            axis.args = list(cex.axis = 3),
            yaxt = "n",
            xaxt = "n",
            xlab = "",
            ylab = ""
        )
        
        xDT<-data.table(xAxAt,xAxLab)
        yDT<-data.table(yAxAt,yAxLab)
        
        axis(1, at = xDT$xAxAt, labels = paste(xDT$xAxLab,"Mb",sep = ""), outer = FALSE, line = NA,lwd.ticks = 2, cex.axis = 2, las = 1)
        #axis(2, at = yDT[yAxLab!=30]$yAxAt, labels = paste(yDT[yAxLab!=30]$yAxLab,"Mb",sep = ""), outer = FALSE, line = NA,lwd.ticks = 2, cex.axis = 2, las = 1)
        
        if (drawGrid) {abline(h = yGrid, v = xGrid, lwd = 2, col = "white")}
        return(matrixToPlot)
    }