## script to look at how entropy scales with area

library(socorro)
library(plyr)

setwd('~/Dropbox/Research/sadScaling')

source('scaling_helpers.R')

paso <- read.csv('../data/stri/PASO.csv', as.is = TRUE)
paso <- paso[paso$year == max(paso$year), ]

## function to calculate singletons (obs and theory)
entr <- function(x, newx) {
    abund <- tapply(newx$count, newx$spp, sum)
    pmf <- as.numeric(table(abund)) / length(abund)
    return(c(S = -sum(pmf*log(pmf))))
}

entropyScale <- scaleMetric(paso, entr)

entropyScale <- ddply(entropyScale, c('scale', 'type'), function(x) {
    c(S = mean(x$S), S.ci = quantile(x$S, c(0.025, 0.975), names = FALSE))
})

pdf('sfiTalk_2017-07/fig_entropyScale.pdf', width = 5, height = 5)
par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2, 0.5, 0),
    bg = 'black', fg = 'white', col.lab = 'white', col.axis = 'white')
with(entropyScale[entropyScale$type == 'spatial', ], {
    plot(scale, S, log = 'xy', axes = FALSE, frame.plot = TRUE, 
         panel.first = segments(x0 = scale, 
                                y0 = ifelse(S.ci1 < 0.01, 0.01, S.ci1), 
                                y1 = S.ci2), 
         ylim = c(0.1, max(S.ci1)), 
         pch = 21, col = 'white', bg = 'black', cex = 1.5, cex.lab = 1.2, 
         xlab = expression('Area ('*m^2*')'))
    logAxis(1, TRUE); logAxis(2, TRUE)
})
dev.off()
