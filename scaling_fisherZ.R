setwd('~/Dropbox/Research/maxentConcept')

library(meteR)
library(parallel)
library(socorro)
library(plyr)

## source hidden function from meteR to deal with areas
source('~/Dropbox/Research/meteR/R/sar_helper_funs.R')

## wd storing the data
dataWD <- '~/Dropbox/Research/data/stri'

## z-values for SADs across scales, for all plots
scaleZ <- mclapply(list.files(dataWD), mc.cores = 6, function(f) {
    x <- read.csv(file.path(dataWD, f), as.is = TRUE)
    ## only look at most recent census
    x <- x[x$year == max(x$year), ]
    
    ## logrithmic intervals for scaling back area
    xscale <- round(max(x$x)) * 2^(-4:0)
    yscale <- round(max(x$y)) * 2^(-4:0)
    
    ## loop over scales and calculate sad
    out <- lapply(1:length(xscale), function(i) {
        newx <- x[x$x <= xscale[i] & x$y <= yscale[i], ]
        thisSAD <- sad(meteESF(newx$spp, newx$count, newx$dbh^2))
        thisZ <- logLikZ(thisSAD, nrep = 499)$z
        
        return(list(sad = thisSAD, z = thisZ))
    })
    
    areas <- xscale * yscale
    names(out) <- paste('area', areas, sep = '')
    
    return(out)
})

names(scaleZ) <- gsub('.csv', '', list.files(dataWD))


## plotting

pdf('fig_scalingSAD.pdf', width = 7, height = 15)
mat <- matrix(as.vector(outer(c(1:length(scaleZ[[1]]), 
                                rep(length(scaleZ[[1]]) + 1, length(scaleZ[[1]]))), 
                              (0:(length(scaleZ) - 1)) * (length(scaleZ[[1]])+1), '+')), 
              nrow = length(scaleZ)*2, byrow = TRUE)

layout(mat, heights = rep(c(2, 1.5), length(scaleZ)))
par(oma = c(2, 0, 1, 2))

for(i in 1:length(scaleZ)) {
    par(mar = rep(0.1, 4))
    for(j in 1:length(scaleZ[[i]])) {
        plot(scaleZ[[i]][[j]]$sad, ptype = 'rad', log = 'y', add.legend = FALSE, 
             axes = FALSE, frame.plot = TRUE)
    }
    
    mtext(names(scaleZ)[i], side = 4, line = 1, xpd = NA)
    
    par(mar = c(3, 0.1, 0.1, 0.1))
    zz <- sapply(scaleZ[[i]], function(x) x$z)
    aa <- as.numeric(gsub('area', '', names(zz)))
    plot(aa, zz, type = 'b', log = 'x', ylim = c(0, ifelse(max(zz) > 5, max(zz), 5)))
    abline(h = 4)
}

dev.off()
