library(raster)
library(MASS)
library(socorro)
library(spectralGP)
library(spatstat)
library(sp)


n <- 2^11

thisGP <- gp(c(n, n), specdens.param = c(0.1, 0.1))
simulate(thisGP)
thisGP <- as.im(0.05 * exp(predict(thisGP)), W = owin(c(0, n), c(0, n)))

thisSSAD <- rpoispp(thisGP)
thisSSAD <- SpatialPoints(cbind(thisSSAD$x, thisSSAD$y))

r <- raster(nrows = n, ncols = n, xmn = 0, xmx = n, ymn = 0, ymx = n)
thisSSAD <- rasterize(thisSSAD, r, fun = 'count')
thisSSAD[is.na(thisSSAD[])] <- 0


negbPar <- matrix(NA, nrow = log(n, base = 2) - 2, ncol = 3)
colnames(negbPar) <- c('A', 'k', 'mu')

for(i in 1:nrow(negbPar)) {
    thisSSAD <- aggregate(thisSSAD, fact = 2, fun = sum)
    f <- fitdistr(thisSSAD[], 'negative binomial')
    
    negbPar[i, 2:3] <- f$estimate
    negbPar[i, 1] <- xmax(extent(thisSSAD)) * ymax(extent(thisSSAD)) / (ncol(thisSSAD) * nrow(thisSSAD))
    
}

plot(negbPar[, c(1, 2)], log = 'xy')
