singPlot <- function(summ, x, y, ...) {
    plot(summ[[x]], summ[[y]], pch = 16,
         col = col, panel.first = {
             lines(summ[[x]], summ[[y]], lty = 2, col = 'gray')
             segments(x0 = summ[[x]], 
                      y0 = summ[[paste(y, 'ci1', sep = '.')]], y1 = summ[[paste(y, 'ci2', sep = '.')]], 
                      col = col)
             segments(x0 = summ[[paste(x, 'ci1', sep = '.')]], x1 = summ[[paste(x, 'ci2', sep = '.')]], 
                      y0 = summ[[y]], 
                      col = col)
         }, 
         xlim = range(summ[, grepl(gsub('perm.', '', x), names(summ))]), 
         ylim = range(summ[, grepl(gsub('perm.', '', y), names(summ))]), 
         ...)
    abline(0, 1)
}

singScaleSumm <- read.csv('scaling_fisherSingletons.csv', as.is = TRUE)

par(mfcol = c(2, 2), mar = c(0, 1, 1, 0) + 0.5, oma = c(3, 2, 0.5, 2), mgp = c(2, 0.25, 0), tcl = -0.25)
singPlot(singScaleSumm[singScaleSumm$site == 'BCI', ], 'n1fish', 'n1obs')
mtext('BCI', side = 3, line = 0.5)
singPlot(singScaleSumm[singScaleSumm$site == 'BCI', ], 'perm.n1fish', 'perm.n1obs', yaxt = 's')

singPlot(singScaleSumm[singScaleSumm$site == 'PASO', ], 'n1fish', 'n1obs')
mtext('PASO', side = 3, line = 0.5)
mtext('Spatial subset', side = 4, outer = FALSE, line = 0.5)
singPlot(singScaleSumm[singScaleSumm$site == 'PASO', ], 'perm.n1fish', 'perm.n1obs', yaxt = 's')
mtext('Random subset', side = 4, outer = FALSE, line = 0.5)

mtext('Logseries-predicted singletons', side = 1, outer = TRUE, line = 1.5)
mtext('Observed singletons', side = 2, outer = TRUE, line = 0.5)





par(mfrow = c(1, 2), mar = c(0, 1, 1, 0) + 0.5, oma = c(3, 2, 0.5, 2), mgp = c(2, 0.25, 0), tcl = -0.25)

singPlot(singScaleSumm[singScaleSumm$site == 'BCI', ], 'n1obs', 'perm.n1obs')
singPlot(singScaleSumm[singScaleSumm$site == 'PASO', ], 'n1obs', 'perm.n1obs')

mtext('Singletons in spatial subsets', side = 1, outer = TRUE, line = 1.5)
mtext('Singletons in random subsets', side = 2, outer = TRUE, line = 0.5)
