library(socorro)

setwd('~/Dropbox/Research/sadScaling')

allZ <- read.csv('scaling_fisherZ.csv', as.is = TRUE)

## function to graph one plot's worth of scaling
plotScaleZ <- function(z, ...) {
    jitterScale <- z$scale*1.2
    
    lwd <- 2
    cex <- 1.5
    permCol <- 'gray60'
    
    plot(z$scale, z$zMean,
         panel.first = {
             lines(jitterScale, z$zPermMean, col = permCol, lwd = lwd)
             segments(x0 = jitterScale, y0 = z$zPermCI1, y1 = z$zPermCI2, col = permCol)
             points(jitterScale, z$zPermMean, 
                    pch = 21, col = permCol, bg = 'white', cex = cex)
             lines(z$scale, z$zMean, lwd = lwd)
             segments(x0 = z$scale, y0 = z$zCI1, y1 = z$zCI2)
         },
         pch = 16, cex = cex, ...)
    
    abline(h = qchisq(0.095, 1), lty = 2, col = 'gray')
}


sortedPlots <- names(sort(tapply(allZ$scale, allZ$plot, max), decreasing = TRUE))
par(mfcol = c(length(sortedPlots), 1), mar = c(0, 0, 0, 2) + 0.5, 
    oma = c(3.5, 3.5, 0, 0), mgp = c(2, 0.75, 0))
for(i in 1:length(sortedPlots)) {
    plotScaleZ(allZ[allZ$plot == sortedPlots[i],], 
               xaxt = 'n', yaxt = 'n', log = 'xy',
               xlim = range(allZ$scale), ylim = c(1e-05, 500), 
               xlab = '', ylab = '')
    axis(2, at = 10^pretty(log(c(1e-05, 500), 10)), 
         labels = sapply(pretty(log(c(1e-05, 500), 10)), 
                         function(p) eval(substitute(expression(10^p), list(p = p)))))
    if(i == length(sortedPlots)) logAxis(1, expLab = TRUE)
    mtext(sortedPlots[i], side = 4, line = 1)
}

mtext(expression('Scale'~(m^2)), side = 1, outer = TRUE, line = 2)
mtext(expression(z^2-'value'), side = 2, outer = TRUE, line = 2)
