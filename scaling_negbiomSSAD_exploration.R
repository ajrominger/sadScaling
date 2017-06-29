setwd('~/Dropbox/Research/sadScaling')

negbScale <- read.csv('scaling_negbiomSSAD.csv', as.is = TRUE)


## exploring how the negative dinomial scales with area. findings so far:
## z-val function is broken
## subsets are poisson (NB fit has large k, so equivilant to pois)
## need to use full fitdistr for NB fitting


getNegB <- function(site, scale, sp) {
    x <- get(tolower(site))
    
    ## make raster object for plot
    r <- raster(ncols = ceiling(max(x$x)), nrows = ceiling(max(x$y)), 
                xmn = 0, xmx = ceiling(max(x$x)), 
                ymn = 0, ymx = ceiling(max(x$y)))
    
    ## determine range of scales
    scales <- determineScale(r)
    
    ## make sure there are minimum 16 cells
    scales <- scales[-(1:4), ]
    
    s <- scales[scales[, 1] * scales[, 2] == scale, ]
    
    res(r) <- s
    
    ## make community matrices for permuted and non-permuted cells
    cells <- cellFromXY(r, xy = x[, c('x', 'y')])
    cellsPerm <- sample(cells) # permuted cells
    
    mat <- tidy2mat(cells, x$spp, x$count)
    permMat <- tidy2mat(cellsPerm, x$spp, x$count)
    
    dat <- mat[, sp]
    permDat <- permMat[, sp]
    
    plot(simpECDF(dat), xlim = range(dat, permDat))
    points(0:max(dat, permDat), 
           pnbinom(0:max(dat, permDat), 
                   size = negbScale$nbPar.k[negbScale$site == site & 
                                                negbScale$spp == sp & 
                                                negbScale$scale == scale], 
                   mu = negbScale$nbPar.mu[negbScale$site == site & 
                                               negbScale$spp == sp & 
                                               negbScale$scale == scale]), 
           col = 'blue', type = 'l')
    points(0:max(dat, permDat), 
           ppois(0:max(dat, permDat), 
                 negbScale$nbPar.mu[negbScale$site == site & 
                                        negbScale$spp == sp & 
                                        negbScale$scale == scale]), 
           col = 'red', type = 'l')
    
    points(simpECDF(permDat), col = 'gray')
    points(0:max(dat, permDat), 
           pnbinom(0:max(dat, permDat), 
                   size = negbScale$perm.nbPar.k[negbScale$site == site & 
                                                     negbScale$spp == sp & 
                                                     negbScale$scale == scale], 
                   mu = negbScale$perm.nbPar.mu[negbScale$site == site & 
                                                    negbScale$spp == sp & 
                                                    negbScale$scale == scale]), 
           col = 'blue', type = 'l', lty = 2)
    points(0:max(dat, permDat), 
           ppois(0:max(dat, permDat), 
                 negbScale$perm.nbPar.mu[negbScale$site == site & 
                                             negbScale$spp == sp & 
                                             negbScale$scale == scale]), 
           col = 'red', type = 'l', lty = 2)
    
    legend('bottomright', legend = paste0('z = ', 
                                          round(negbScale$nbZ[negbScale$site == site & 
                                                                  negbScale$spp == sp & 
                                                                  negbScale$scale == scale], 3)))
}









layout(matrix(c(1, 2, 3, 4, 5, 5), ncol = 3))
par(mar = c(2.5, 2.5, 0, 0) + 0.2, mgp = c(1.5, 0.5, 0))

j <- 1
with(negbScale[negbScale$site == 'PASO', ], {
    i <- unique(spp)[j]
    browser()
    plot(1, xlim = range(scale), ylim = c(min(negbScale$nbLL - negbScale$pLL, na.rm = TRUE), 400), 
         log = 'x', type = 'n', 
         xlab = 'scale', ylab = 'NB logLik - Pois logLik')
    
    points(scale[spp == i], nbLL[spp == i] - pLL[spp == i], type = 'l')
    points(scale[spp == i], perm.nbLL[spp == i] - perm.pLL[spp == i], type = 'l', col = 'red')
    
    plot(scale[spp == i], nbZ[spp == i], log = 'xy', type = 'b', ylim = c(0.0001, 5000))
    points(scale[spp == i], perm.nbZ[spp == i], col = 'red', type = 'b')
    abline(h = qchisq(0.95, 1))
    
    plot(1, xlim = range(scale), ylim = c(0.001, 100), log = 'xy', type = 'n', 
         xlab = 'scale', ylab = 'NB k')
    points(scale[spp == i], nbPar.k[spp == i], type = 'l')
    
    plot(1, xlim = range(scale), ylim = c(0.001, 100), log = 'xy', type = 'n', 
         xlab = 'scale', ylab = 'NB mu')
    points(scale[spp == i], nbPar.mu[spp == i], type = 'l')
    
    plot(1, xlim = range(scale), ylim = c(0.001, 20), log = 'x', type = 'n', 
         xlab = 'scale', ylab = 'KL')
    points(scale[spp == i], kl[spp == i], type = 'l')
    points(scale[spp == i], perm.kl[spp == i], type = 'l', col = 'red')
    
    # getNegB('PASO', unique(negbScale$scale)[5], i)
})
j <- j + 1

layout(matrix(1:12, nrow = 3))
par(mar = c(2, 2, 0, 0))
lapply(unique(negbScale$scale), function(foo) {
    getNegB('PASO', foo, i)
})
