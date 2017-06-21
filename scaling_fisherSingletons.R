library(socorro)
library(plyr)
singleton <- function(x, newx) {
    abund <- tapply(newx$count, newx$spp, sum)
    n1fish <- sum(sad2Rank(sad(abund, model = 'fish')) == 1)
    n1obs <- sum(abund == 1)
    return(c(n1fish = n1fish, n1obs = n1obs))
}

foo <- scaleMetric(bci, singleton)
bla <- ddply(foo, 'scale', function(x) {
    out <- c(colMeans(x[, -1]), 
             as.vector(apply(x[, -1], 2, quantile, probs = c(0.025, 0.975))))
    names(out) <- c(names(out)[1:4], 
                    paste(rep(names(out)[1:4], each = 2), c('ci1', 'ci2'), sep = '.'))
    return(out)
})

par(mfrow = c(1, 2), mar = c(0, 0, 2, 0) + 0.5, oma = c(3, 3, 1, 1), mgp = c(2, 0.75, 0))
plot(bla$n1fish, bla$n1obs, pch = 16,
     col = quantCol(bla$scale, pal = hsv(c(0.12, 0.6, 1), c(0.8, 1, 1), c(0.8, 0.8, 0.7)), 
                    trans = 'log'), 
     panel.first = {
         lines(bla$n1fish, bla$n1obs, lty = 2, col = 'gray')
         # segments(x0 = bla$n1fish, y0 = bla$n1obs.ci1, y1 = bla$n1obs.ci2, 
         #          col = quantCol(bla$scale, pal = hsv(c(0.2, 0.6, 1)), 
         #                         trans = 'log'))
         # segments(x0 = bla$n1fish.ci1, x1 = bla$n1fish.ci2, y0 = bla$n1obs, 
         #          col = quantCol(bla$scale, pal = c('yellow', 'blue', 'red'), 
         #                         trans = 'log'))
     }, 
     xlim = range(bla[, grepl('n1fish', names(bla))]), 
     ylim = range(bla[, grepl('n1obs', names(bla))]))
abline(0, 1)
mtext('Spatial subset', side = 3, outer = FALSE, line = 0)

plot(bla$perm.n1fish, bla$perm.n1obs, pch = 21, bg = 'white',
     col = quantCol(bla$scale, pal = c('yellow', 'blue', 'red'), 
                    trans = 'log'), 
     panel.first = {
         lines(bla$perm.n1fish, bla$perm.n1obs, lty = 2, col = 'gray')
         segments(x0 = bla$perm.n1fish, y0 = bla$perm.n1obs.ci1, y1 = bla$perm.n1obs.ci2, 
                  col = quantCol(bla$scale, pal = c('yellow', 'blue', 'red'), 
                                 trans = 'log'))
         segments(x0 = bla$perm.n1fish.ci1, x1 = bla$perm.n1fish.ci2, y0 = bla$perm.n1obs, 
                  col = quantCol(bla$scale, pal = c('yellow', 'blue', 'red'), 
                                 trans = 'log'))
     }, 
     yaxt = 'n',
     xlim = range(bla[, grepl('perm.n1fish', names(bla))]), 
     ylim = range(bla[, grepl('perm.n1obs', names(bla))]))
abline(0, 1)

mtext('Random subset', side = 3, outer = FALSE, line = 0)

mtext('Logseries-predicted singletons', side = 1, outer = TRUE, line = 1.5)
mtext('Observed singletons', side = 2, outer = TRUE, line = 1.5)
