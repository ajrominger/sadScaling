library(raster)
library(MASS)
library(socorro)
library(spectralGP)
library(spatstat)
library(sp)


# make a spatial GP on a fine grid
n <- 2^8

thisGP <- gp(c(n, n), specdens.param = c(0.1, 1))
simulate(thisGP)
plot(as.im(predict(thisGP)))


relGP <- as.im(0.001 * exp(predict(thisGP)), W = owin(c(0, n), c(0, n)))

# sample it as a Poisson
thisSSAD <- rpoispp(relGP)
thisSSAD <- SpatialPoints(cbind(thisSSAD$x, thisSSAD$y))

nrow(thisSSAD@coords)

pdf('foo.pdf')
plot(log(relGP))
points(thisSSAD, pch = 16, col = hsv(0, 0, 1, alpha = 1), cex = 0.2)
dev.off()


# make a raster of counts from the Poisson process
r <- raster(nrows = n / 2, ncols = n / 2, xmn = 0, xmx = n, ymn = 0, ymx = n)
thisSSAD <- rasterize(thisSSAD, r, fun = 'count', background = 0)


# calculate negbinom pars across increasing aggregations
negbPar <- matrix(NA, nrow = log(n, base = 2) - 2, ncol = 3)
colnames(negbPar) <- c('A', 'k', 'mu')


for(i in 1:nrow(negbPar)) {
    if(i > 1) {
        thisSSAD <- aggregate(thisSSAD, fact = 2, fun = sum)
    }
    
    f <- fitdistr(thisSSAD[], 'negative binomial')
    
    negbPar[i, 2:3] <- f$estimate
    negbPar[i, 1] <- xmax(extent(thisSSAD)) * ymax(extent(thisSSAD)) / 
        (ncol(thisSSAD) * nrow(thisSSAD))
    
}

plot(negbPar[, c(1, 2)], log = 'xy', axes = FALSE, frame.plot = TRUE, ylim = c(1, 10))
logAxis(1, expLab = TRUE)
logAxis(2)
