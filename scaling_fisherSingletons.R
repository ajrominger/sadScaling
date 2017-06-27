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
singleton <- function(x, newx) {
    abund <- tapply(newx$count, newx$spp, sum)
    n1fish <- sum(sad2Rank(sad(abund, model = 'fish')) == 1)
    n1obs <- sum(abund == 1)
    return(c(n1fish = n1fish, n1obs = n1obs))
}

## calculate singletons across scales
singScaleBCI <- scaleMetric(bci, singleton)
singScalePASO <- scaleMetric(paso, singleton)
singScale <- rbind(cbind(site = 'BCI', singScaleBCI), cbind(site = 'PASO', singScalePASO))

singScaleSumm <- ddply(singScale, c('site', 'scale'), function(x) {
    out <- c(colMeans(x[, -(1:2)]), 
             n1diff = mean(x[, 3] - x[, 4]),
             perm.n1diff = mean(x[, 5] - x[, 6]),
             as.vector(apply(x[, -(1:2)], 2, quantile, probs = c(0.025, 0.975))), 
             quantile(x[, 3] - x[, 4], probs = c(0.025, 0.975)), 
             quantile(x[, 5] - x[, 6], probs = c(0.025, 0.975)))
    names(out) <- c(names(out)[1:6], 
                    paste(rep(names(out)[1:6], each = 2), c('ci1', 'ci2'), sep = '.'))
    return(out)
})

## write out
write.csv(singScaleSumm, 'scaling_fisherSingletons.csv', row.names = FALSE)
