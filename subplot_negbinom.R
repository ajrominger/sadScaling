library(MASS)
library(plyr)
library(socorro)

## source hidden function from meteR to deal with areas
source('~/Dropbox/Research/meteR/R/sar_helper_funs.R')

## wd storing the data
dataWD <- '~/Dropbox/Research/data/stri'

## helper function to make all calculations
ssadLL <- function(f, nrow, ncol) {
    # browser()
    dat <- read.csv(file.path(dataWD, f))
    dat <- dat[dat$year == max(dat$year), ]
    datArea <- .findAreas(dat$spp, dat$count, x = dat$x, y = dat$y, row = nrow, col = ncol)
    dat$cell <- paste(datArea$row, datArea$col, sep = ',')
    
    datMat <- tidy2mat(dat$cell, dat$spp, dat$count)
    datMat <- datMat[, colSums(datMat) > 0]
    
    # datMat2 <- datMat[round(seq(1, 50, length.out = 20)), ]
    # datMat2 <- datMat2[, colSums(datMat2) > 0]
    # 
    # llObs2 <- apply(datMat2, 2, function(x) {
    #     unlist(fitdistr(x, 'negative binomial')[c('estimate', 'loglik')])
    # })
    
    llObs <- apply(datMat, 2, function(x) {
        unlist(fitdistr(x, 'negative binomial')[c('estimate', 'loglik')])
    })
    
    # foo2 <- llObs2[2, ]
    # foo <- llObs[2, names(foo2)]
    # 
    # plot(foo, foo2)
    
    llThr <- apply(llObs, 2, function(x) 
        sum(dnbinom(rnbinom(100*nrow(datMat), size = x[1], mu = x[2]), 
                    size = x[1], mu = x[2], log = TRUE)) / 100)
    
    return(data.frame(t(llObs), llThr = llThr, N = colSums(datMat)))
}


## BCI ssad...all negative binomial
bciSSAD <- ssadLL('BCIS.csv', 5, 10)
plot(bciSSAD[, 3:4])
abline(0, 1, col = 'red')

## UCSC...all negative binomial
ucscSSAD <- ssadLL('UCSC.csv', 6, 4)
plot(ucscSSAD[, 3:4])
abline(0, 1, col = 'red')

## Cocoli...all negative binomial
cocoSSAD <- ssadLL('COCO.csv', 6, 2)
plot(cocoSSAD[, 3:4])
abline(0, 1, col = 'red')

## Paso...all negative binomial
pasoSSAD <- ssadLL('PASO.csv', 5, 10)
plot(pasoSSAD[, 3:4])
abline(0, 1, col = 'red')


plot(sort(pasoSSAD[, 1], TRUE), type = 'l', log = 'y', 
     ylim = range(pasoSSAD[, 1], bciSSAD[, 1]), col = 'red', lwd = 2)
points(sort(bciSSAD[, 1], TRUE), type = 'l', col = 'blue', lwd = 2)
points(sort(cocoSSAD[, 1], TRUE), type = 'l', col = 'skyblue', lwd = 2)
points(sort(ucscSSAD[, 1], TRUE), type = 'l', col = 'black', lwd = 2)

plot(sort(pasoSSAD[, 2], TRUE), type = 'l', log = 'y', 
     ylim = range(pasoSSAD[, 2], bciSSAD[, 2]), col = 'red', lwd = 2)
points(sort(bciSSAD[, 2], TRUE), type = 'l', col = 'blue', lwd = 2)
points(sort(cocoSSAD[, 2], TRUE), type = 'l', col = 'skyblue', lwd = 2)
points(sort(ucscSSAD[, 2], TRUE), type = 'l', col = 'black', lwd = 2)

plot(pasoSSAD$estimate.mu, pasoSSAD$N)



pasoSim <- replicate(50, rnbinom(nrow(pasoSSAD), size = pasoSSAD$estimate.size, 
                                 mu = pasoSSAD$estimate.mu))

plot(sad(meteESF(1:nrow(pasoSim), rowSums(pasoSim[, 3, drop = FALSE]))), 
     ptype = 'rad', log = 'y')
