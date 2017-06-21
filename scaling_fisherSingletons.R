## script to look at how singletons scale with area and how 
## well the log series predicts them

library(socorro)
library(plyr)

setwd('~/Dropbox/Research/sadScaling')

source('scaling_helpers.R')

bci <- read.csv('../data/stri/BCIS.csv', as.is = TRUE)
bci <- bci[bci$year == max(bci$year), ]

## function to calculate singletons (obs and theory)
singleton <- function(x, newx) {
    abund <- tapply(newx$count, newx$spp, sum)
    n1fish <- sum(sad2Rank(sad(abund, model = 'fish')) == 1)
    n1obs <- sum(abund == 1)
    return(c(n1fish = n1fish, n1obs = n1obs))
}

## calculate singletons across scales
singScale <- scaleMetric(bci, singleton)

singScaleSumm <- ddply(singScale, 'scale', function(x) {
    out <- c(colMeans(x[, -1]), 
             as.vector(apply(x[, -1], 2, quantile, probs = c(0.025, 0.975))))
    names(out) <- c(names(out)[1:4], 
                    paste(rep(names(out)[1:4], each = 2), c('ci1', 'ci2'), sep = '.'))
    return(out)
})

par(mfrow = c(1, 2), mar = c(0, 0, 2, 0) + 0.5, 
    oma = c(3, 3, 1, 1), mgp = c(2, 0.75, 0))
col <- quantCol(singScaleSumm$scale, pal = hsv(c(0.12, 0.6, 1), c(0.8, 1, 1), 
                                               c(0.8, 0.8, 0.7)), 
                trans = 'log')
plot(singScaleSumm$n1fish, singScaleSumm$n1obs, pch = 16,
     col = col, panel.first = {
         lines(singScaleSumm$n1fish, singScaleSumm$n1obs, lty = 2, col = 'gray')
         segments(x0 = singScaleSumm$n1fish, 
                  y0 = singScaleSumm$n1obs.ci1, y1 = singScaleSumm$n1obs.ci2, 
                  col = col)
         segments(x0 = singScaleSumm$n1fish.ci1, x1 = singScaleSumm$n1fish.ci2, 
                  y0 = singScaleSumm$n1obs, 
                  col = col)
     }, 
     xlim = range(singScaleSumm[, grepl('n1fish', names(singScaleSumm))]), 
     ylim = range(singScaleSumm[, grepl('n1obs', names(singScaleSumm))]))
abline(0, 1)
mtext('Spatial subset', side = 3, outer = FALSE, line = 0)

plot(singScaleSumm$perm.n1fish, singScaleSumm$perm.n1obs, pch = 16,
     col = col, 
     panel.first = {
         lines(singScaleSumm$perm.n1fish, singScaleSumm$perm.n1obs, lty = 2, col = 'gray')
         segments(x0 = singScaleSumm$perm.n1fish, 
                  y0 = singScaleSumm$perm.n1obs.ci1, y1 = singScaleSumm$perm.n1obs.ci2, 
                  col = col)
         segments(x0 = singScaleSumm$perm.n1fish.ci1, x1 = singScaleSumm$perm.n1fish.ci2, 
                  y0 = singScaleSumm$perm.n1obs, 
                  col = col)
     }, 
     yaxt = 'n',
     xlim = range(singScaleSumm[, grepl('perm.n1fish', names(singScaleSumm))]), 
     ylim = range(singScaleSumm[, grepl('perm.n1obs', names(singScaleSumm))]))
abline(0, 1)

mtext('Random subset', side = 3, outer = FALSE, line = 0)

mtext('Logseries-predicted singletons', side = 1, outer = TRUE, line = 1.5)
mtext('Observed singletons', side = 2, outer = TRUE, line = 1.5)
