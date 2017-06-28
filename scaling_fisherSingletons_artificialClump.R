## script to look at how singletons scale with area when species geographic distributions
## are artificially forced to be clumped

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

## make clumped geographic distribuitons
x <- bci
temp <- x$x + ceiling(x$y / 50) * 50
x$x <- x$x[order(temp)]
x$y <- x$y[order(temp)]
x$spp <- x$spp[unlist(lapply(names(sort(tapply(x$count, x$spp, sum))), function(spp) which(x$spp == spp)))]
x$spp <- x$spp[order(tapply(x$count, x$spp, sum)[x$spp])]

singScaleSim <- scaleMetric(x, singleton)

foo <- ddply(singScaleSim, 'scale', function(x) {
    out <- c(colMeans(x[, -1]), 
             n1diff = mean(x[, 2] - x[, 3]),
             perm.n1diff = mean(x[, 4] - x[, 5]),
             as.vector(apply(x[, -1], 2, quantile, probs = c(0.025, 0.975))), 
             quantile(x[, 2] - x[, 3], probs = c(0.025, 0.975)), 
             quantile(x[, 4] - x[, 5], probs = c(0.025, 0.975)))
    names(out) <- c(names(out)[1:6], 
                    paste(rep(names(out)[1:6], each = 2), c('ci1', 'ci2'), sep = '.'))
    return(out)
})


par(mfcol = 1:2, mar = c(2, 2, 0, 0))
scalePlot(singScaleSumm[singScaleSumm$site == 'BCI', ], 'n1obs', log = 'x', ylim = c(0, 50))
scalePlot(singScaleSumm[singScaleSumm$site == 'BCI', ], 'perm.n1obs', log = 'x', col = 'red', add = TRUE)
# abline(h = 0, col = 'gray')

scalePlot(foo, 'n1obs', log = 'x', ylim = c(0, 50))
scalePlot(foo, 'perm.n1obs', log = 'x', col = 'red', add = TRUE)
abline(h = 0, col = 'gray')
