scalePlot <- function(summ, y, jit = 1, col = 'black', ylim = NULL, add = FALSE, ...) {
    x <- summ$scale * jit
    
    if(is.null(ylim)) {
        ylim <- range(summ[, grep(y, colnames(summ))])
    }
    
    if(add) {
        lines(x, summ[[y]], col = col, lty = 2)
        segments(x0 = x, 
                 y0 = summ[[paste(y, 'ci1', sep = '.')]], y1 = summ[[paste(y, 'ci2', sep = '.')]], 
                 col = col)
        points(x, summ[[y]], pch = 16, col = col)
    } else {
        plot(x, summ[[y]], pch = 16, col = col, ylim = ylim,
             panel.first = {
                 lines(x, summ[[y]], col = col, lty = 2)
                 segments(x0 = x, 
                          y0 = summ[[paste(y, 'ci1', sep = '.')]], y1 = summ[[paste(y, 'ci2', sep = '.')]], 
                          col = col)
             }, ...)
    }
}

singScaleSumm <- read.csv('scaling_fisherSingletons.csv', as.is = TRUE)
relLLScaleSumm <- read.csv('scaling_fisherRelLL.csv', as.is = TRUE)
zScaleSumm <- read.csv('scaling_fisherZ.csv', as.is = TRUE)

obsCol <- hsv(0.5, 1, 0.5)
fishCol <- hsv(0, 1, 0.5)
perm.obsCol <- hsv(0.6, 0.9, 0.6)
perm.fishCol <- hsv(0.1, 0.8, 0.7)

par(mfrow = c(2, 2), mar = c(0, 0, 0, 0) + 0.5, oma = c(3, 3, 2, 2), mgp = c(2, 0.25, 0), tcl = -0.25)

scalePlot(singScaleSumm[singScaleSumm$site == 'BCI', ], 'n1obs', log = 'x', xaxt = 'n', col = obsCol)
scalePlot(singScaleSumm[singScaleSumm$site == 'BCI', ], 'n1fish', add = TRUE, jit = 1.1, col = fishCol)
legend('topleft', legend = c('Observed', 'Logseries'), 
       col = c(obsCol, fishCol), pch = 16, bty = 'n')
mtext('Spatial subsets', side = 3, line = 1)

scalePlot(singScaleSumm[singScaleSumm$site == 'BCI', ], 'perm.n1obs', log = 'x', xaxt = 'n', yaxt = 'n', 
          col = perm.obsCol)
scalePlot(singScaleSumm[singScaleSumm$site == 'BCI', ], 'perm.n1fish', add = TRUE, jit = 1.1, 
          col = perm.fishCol)
legend('topleft', legend = c('Observed', 'Logseries'), 
       col = c(perm.obsCol, perm.fishCol), pch = 16, bty = 'n')
mtext('BCI', side = 4, line = 0.5)
mtext('Random subsets', side = 3, line = 1)

scalePlot(singScaleSumm[singScaleSumm$site == 'PASO', ], 'n1obs', log = 'x', xaxt = 'n', col = obsCol)
scalePlot(singScaleSumm[singScaleSumm$site == 'PASO', ], 'n1fish', add = TRUE, jit = 1.1, col = fishCol, 
          xaxt = 'n')
logAxis(1, expLab = TRUE)

scalePlot(singScaleSumm[singScaleSumm$site == 'PASO', ], 'perm.n1obs', log = 'x', xaxt = 'n', yaxt = 'n', 
          col = perm.obsCol)
scalePlot(singScaleSumm[singScaleSumm$site == 'PASO', ], 'perm.n1fish', add = TRUE, jit = 1.1,
          xaxt = 'n', col = perm.fishCol)
logAxis(1, expLab = TRUE)
mtext('PASO', side = 4, line = 0.5)

mtext(expression('Area ('*m^2*')'), side = 1, line = 1, outer = TRUE)
mtext('Number of singletons', side = 2, line = 1, outer = TRUE)


par(mfrow = c(2, 1), mar = c(0, 0, 0, 0) + 0.5, oma = c(3, 3, 2, 2), mgp = c(2, 0.25, 0), tcl = -0.25)

scalePlot(singScaleSumm[singScaleSumm$site == 'BCI', ], 'n1obs', log = 'x', xaxt = 'n', col = obsCol)
scalePlot(singScaleSumm[singScaleSumm$site == 'BCI', ], 'perm.n1obs', add = TRUE, jit = 1.1, col = perm.obsCol)
mtext('BCI', side = 4, line = 0.5)
legend('topleft', legend = c('Spatial', 'Random'), col = c(obsCol, perm.obsCol), pch = 16, bty = 'n')

scalePlot(singScaleSumm[singScaleSumm$site == 'PASO', ], 'n1obs', log = 'x', xaxt = 'n', col = obsCol)
scalePlot(singScaleSumm[singScaleSumm$site == 'PASO', ], 'perm.n1obs', add = TRUE, jit = 1.1, col = perm.obsCol)
logAxis(1, expLab = TRUE)
mtext('PASO', side = 4, line = 0.5)

mtext(expression('Area ('*m^2*')'), side = 1, line = 1, outer = TRUE)
mtext('Observed number of singletons', side = 2, line = 1, outer = TRUE)


par(mfrow = c(2, 1), mar = c(0, 0, 0, 0) + 0.5, oma = c(3, 3, 2, 2), mgp = c(2, 0.25, 0), tcl = -0.25)

scalePlot(zScaleSumm[zScaleSumm$site == 'BCI', ], 'z2', log = 'x', xaxt = 'n', col = fishCol)
scalePlot(zScaleSumm[zScaleSumm$site == 'BCI', ], 'perm.z2', col = perm.fishCol, add = TRUE, jit = 1.1)
abline(h = qchisq(0.95, 1), col = 'gray')

scalePlot(zScaleSumm[zScaleSumm$site == 'PASO', ], 'z2', log = 'x', xaxt = 'n', col = fishCol)
scalePlot(zScaleSumm[zScaleSumm$site == 'PASO', ], 'perm.z2', col = perm.fishCol, add = TRUE, jit = 1.1)
abline(h = qchisq(0.95, 1), col = 'gray')
logAxis(1, expLab = TRUE)

mtext(expression('Area ('*m^2*')'), side = 1, line = 1, outer = TRUE)
mtext(expression(z^2*'-score'), side = 2, outer = TRUE, line = 1)




par(mfrow = c(2, 2), mar = c(0, 0, 0, 0) + 0.5, oma = c(3, 3, 2, 2), mgp = c(2, 0.25, 0), tcl = -0.25)
scalePlot(relLLScaleSumm[relLLScaleSumm$site == 'BCI', ], 'llST', log = 'x')
scalePlot(relLLScaleSumm[relLLScaleSumm$site == 'BCI', ], 'llF', col = 'red', add = TRUE)

scalePlot(relLLScaleSumm[relLLScaleSumm$site == 'BCI', ], 'perm.llST', log = 'x')
scalePlot(relLLScaleSumm[relLLScaleSumm$site == 'BCI', ], 'perm.llF', col = 'red', add = TRUE)

scalePlot(relLLScaleSumm[relLLScaleSumm$site == 'PASO', ], 'llST', log = 'x')
scalePlot(relLLScaleSumm[relLLScaleSumm$site == 'PASO', ], 'llF', col = 'red', add = TRUE)

scalePlot(relLLScaleSumm[relLLScaleSumm$site == 'PASO', ], 'perm.llST', log = 'x')
scalePlot(relLLScaleSumm[relLLScaleSumm$site == 'PASO', ], 'perm.llF', col = 'red', add = TRUE)








par(mfcol = c(3, 1), mar = c(2, 2, 0, 0))

scalePlot(singScaleSumm[singScaleSumm$site == 'PASO', ], 'n1diff', log = 'x', xaxt = 'n')
scalePlot(singScaleSumm[singScaleSumm$site == 'PASO', ], 'perm.n1diff', log = 'x', xaxt = 'n', col = 'red', add = TRUE)
abline(h = 0, col = 'gray')

scalePlot(relLLScaleSumm[relLLScaleSumm$site == 'PASO', ], 'relLL', log = 'x')
scalePlot(relLLScaleSumm[relLLScaleSumm$site == 'PASO', ], 'perm.relLL', log = 'x', col = 'red', add = TRUE)
abline(h = 0, col = 'gray')






par(mfcol = c(3, 1), mar = c(2, 2, 0, 0))
scalePlot(singScaleSumm[singScaleSumm$site == 'BCI', ], 'n1diff', log = 'x', xaxt = 'n')
scalePlot(singScaleSumm[singScaleSumm$site == 'BCI', ], 'perm.n1diff', log = 'x', xaxt = 'n', col = 'red', add = TRUE)
abline(h = 0, col = 'gray')

scalePlot(relLLScaleSumm[relLLScaleSumm$site == 'BCI', ], 'relLL', log = 'x')
scalePlot(relLLScaleSumm[relLLScaleSumm$site == 'BCI', ], 'perm.relLL', log = 'x', col = 'red', add = TRUE)
abline(h = 0, col = 'gray')

scalePlot(zScaleSumm[zScaleSumm$site == 'BCI', ], 'z2', log = 'x')
scalePlot(zScaleSumm[zScaleSumm$site == 'BCI', ], 'perm.z2', log = 'x', col = 'red', add = TRUE)
abline(h = qchisq(0.95, 1), col = 'gray')
