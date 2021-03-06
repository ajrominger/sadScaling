## script to look at how singletons scale with area and how 
## well the log series predicts them

library(socorro)
library(plyr)

setwd('~/Dropbox/Research/sadScaling')

source('scaling_helpers.R')

bci <- read.csv('../data/stri/BCIS.csv', as.is = TRUE)
bci <- bci[bci$year == max(bci$year), ]

paso <- read.csv('../data/stri/PASO.csv', as.is = TRUE)
paso <- paso[paso$year == max(paso$year), ]

## function to calculate singletons (obs and theory)
z2 <- function(x, newx) {
    abund <- tapply(newx$count, newx$spp, sum)
    thisSAD <- sad(abund, model = 'fish', keepData = TRUE)
    return(c(z2 = logLikZ(thisSAD)$z))
}

## calculate singletons across scales
zScaleBCI <- scaleMetric(bci, z2)
zScalePASO <- scaleMetric(paso, z2)
zScale <- rbind(cbind(site = 'BCI', zScaleBCI), cbind(site = 'PASO', zScalePASO))

zScaleSumm <- ddply(zScale, c('site', 'scale', 'type'), function(x) {
    out <- c(colMeans(x[, c('N', 'z2')]), 
             z2.ci = quantile(x$z2, probs = c(0.025, 0.975), names = FALSE))
    return(out)
})

## write out
write.csv(zScaleSumm, 'scaling_fisherZ_NEW.csv', row.names = FALSE)
