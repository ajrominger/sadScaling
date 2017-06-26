## script to look at how logseries logLik and sample transform logLik
## scale with area

library(socorro)
library(plyr)

setwd('~/Dropbox/Research/sadScaling')

source('../happySAD/sampleTransform.R')

bci <- read.csv('../data/stri/BCIS.csv', as.is = TRUE)
bci <- bci[bci$year == max(bci$year), ]

paso <- read.csv('../data/stri/PASO.csv', as.is = TRUE)
paso <- paso[paso$year == max(paso$year), ]

library(pika)
library(raster)
library(parallel)

## function to take plot data, divide it up across scales, and caluclate
## metrics across scales
#' @param x is the data
#' @param fun is a function returning the metric of interest
#' @param perm is a logical indicating whether to compute the metric on permuted cells

scaleMetric <- function(x, fun, perm = TRUE) {
    ## determine which dimension to scale
    if(max(x$x) >= max(x$y)) {
        scaleVar <- 'x'
    } else {
        scaleVar <- 'y'
    }
    
    ## make sequence of scales along which to cut
    varMax <- ceiling(max(x[[scaleVar]]))
    scales <- seq(varMax, varMax/2, length.out = 6)
    
    ## loop over scales, cutting data by re-scaled raster and calculating SAD z-values
    out <- mclapply(1:length(scales), mc.cores = 6, FUN = function(s) {
        ## get subsetted data
        theseInd <- x[[scaleVar]] <= scales[s]
        theseIndPerm <- sample(theseInd)
        
        m <- calcAtCell(x, theseInd, fun)
        
        if(s == 1 | !perm) {
            mPerm <- m
        } else {
            mPerm <- calcAtCell(x, theseIndPerm, fun)
        }
        
        return(c(m, perm = mPerm))
    })
        
    allM <- do.call(rbind, out)
    
    ## return scale and metric values
    return(data.frame(scale = scales * ifelse(scaleVar == 'x', ceiling(max(x$y)), ceiling(max(x$x))), allM))
}

## helper function to subset data by cell and calculate metric within that cell
#' @param x is the data
#' @param inds individuals to sample
#' @param fun is the function to calculate the metric
calcAtCell <- function(x, inds, fun) {
    newx <- x[inds, ]
    return(fun(x, newx))
}

## function to calculate log lik under fisher and under sample transform
relLL <- function(x, newx) {
    abund <- tapply(newx$count, newx$spp, sum)
    llF <- logLik(sad(abund, 'fish'))
    
    temp <- sampTrans(abund, sad(tapply(x$count, x$spp, sum)), nrow(newx), 
                      include0 = FALSE, log = TRUE)
    temp[!is.finite(temp)] <- min(temp[is.finite(temp)])
    llST <- sum(temp)
    
    return(c(llF = llF, llST = llST))
}

## calculate singletons across scales
relLLScaleBCI <- scaleMetric(bci, relLL)
relLLScalePASO <- scaleMetric(paso, relLL)
relLLScale <- rbind(cbind(site = 'BCI', relLLScaleBCI), cbind(site = 'PASO', relLLScalePASO))
