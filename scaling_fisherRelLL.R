## script to look at how logseries logLik and sample transform logLik
## scale with area

library(socorro)
library(plyr)

setwd('~/Dropbox/Research/sadScaling')

source('scaling_helpers.R')
source('../happySAD/sampleTransform.R')

bci <- read.csv('../data/stri/BCIS.csv', as.is = TRUE)
bci <- bci[bci$year == max(bci$year), ]

paso <- read.csv('../data/stri/PASO.csv', as.is = TRUE)
paso <- paso[paso$year == max(paso$year), ]

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

relLLScaleSumm <- ddply(relLLScale, c('site', 'scale'), function(x) {
    out <- c(colMeans(x[, -(1:2)]), 
             relLL = mean(x[, 3] - x[, 4]), 
             perm.relLL = mean(x[, 5] - x[, 6]),
             as.vector(apply(x[, -(1:2)], 2, quantile, probs = c(0.025, 0.975))), 
             quantile(x[, 3] - x[, 4], probs = c(0.25, 0.975)), 
             quantile(x[, 5] - x[, 6], probs = c(0.25, 0.975)))
    names(out) <- c(names(out)[1:6], 
                    paste(rep(names(out)[1:6], each = 2), c('ci1', 'ci2'), sep = '.'))
    return(out)
})

## write out
write.csv(relLLScaleSumm, 'scaling_fisherRelLL.csv', row.names = FALSE)
