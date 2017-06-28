## script to look at how singletons scale with area and how 
## well the log series predicts them

library(socorro)
library(plyr)
library(MASS)
library(pika)
library(raster)
library(parallel)

setwd('~/Dropbox/Research/sadScaling')

source('scaling_helpers.R')

bci <- read.csv('../data/stri/BCIS.csv', as.is = TRUE)
bci <- bci[bci$year == max(bci$year), ]

paso <- read.csv('../data/stri/PASO.csv', as.is = TRUE)
paso <- paso[paso$year == max(paso$year), ]

foo <- replicate(500, {
    bci1 <- bci[sample(nrow(bci), 10000), ]
    sum(bci1$spp == 'CROTBI')
})

foofit <- fitdistr(foo, 'negative binomial')
plot(simpECDF(foo))
curve(pnbinom(x, size = foofit$estimate[1], mu = foofit$estimate[2]), add = TRUE, col = 'red')
curve(ppois(x, mean(foo)), add = TRUE, col = 'blue')

## fit negbinom and pois and compare logLik of the two
## figure out sampTrans type thing for neg binom
## compare KL distance of best fit nebginom and pois



## function to take plot data, divide it up across scales, and caluclate
## metrics across scales
#' @param x is the data
#' @param fun is a function returning the metric of interest
#' @param perm is a logical indicating whether to compute the metric on permuted cells

scaleNegBinom <- function(x, perm = TRUE) {
    ## make raster object for plot
    r <- raster(ncols = ceiling(max(x$x)), nrows = ceiling(max(x$y)), 
                xmn = 0, xmx = ceiling(max(x$x)), 
                ymn = 0, ymx = ceiling(max(x$y)))
    
    ## determine range of scales
    scales <- determineScale(r)
    
    ## make sure there are minimum 16 cells
    scales <- scales[-(1:4), ]
    
    ## loop over scales, cutting data by re-scaled raster and calculating SAD z-values
    out <- lapply(1:nrow(scales), FUN = function(s) {
        ## rescale
        res(r) <- scales[s, ]
        
        ## make community matrices for permuted and non-permuted cells
        cells <- cellFromXY(r, xy = x[, c('x', 'y')])
        cellsPerm <- sample(cells) # permuted cells
        mat <- tidy2mat(cells, x$spp, x$count)
        if(perm) permMat <- tidy2mat(cells, x$spp, x$count)
        
        sppOut <- lapply(1:ncol(mat), function(i) {
            ## get distrib params
            mu <- mean(foo)
            size <- mu^2 / (var(foo) - mu)
            
            nbPar <- c(k = mean(mat[, i]) / (var(mat[, i]) - mean(mat[, i])), 
                       mu = mean(mat[, i]))
            pPar <- mean(mat[, i])
            
            ## calculate likelihoods
            nbLL <- sum(dnbinom(mat[, i], nbPar[1], mu = nbPar[2], log = TRUE))
            pLL <- sum(dpois(mat[, i], pPar, log = TRUE))
            
            ## calculate z2 for negbinom
            p0 <- dnbinom(0:10^6, size = nbPar[1], mu = nbPar[2], log = TRUE)
            p0 <- p0[is.finite(p0)]
            n <- nrow(mat)
            m <- sum(p0 * exp(p0)) * n
            v <- sum((m/n - p0)^2 * exp(p0)) * n
            nbZ <- ((nbLL - m) / sqrt(v))^2
            
            ## calculate KL dist for neg binom and pois
            q0 <- dpois(seq(0, length.out = length(p0)), pPar, log = TRUE)
            p0 <- p0[is.finite(q0)]
            q0 <- q0[is.finite(q0)]
            kl <- sum(exp(p0) * (p0 - q0))
            
            return(c(nbLL = nbLL, pLL = pLL, nbZ = nbZ, nbPar = nbPar, kl = kl))
        })
        
        sppOut <- do.call(rbind, sppOut)
        
        ## return scale and metric values
        return(data.frame(scale = prod(scales[s, ]), spp = colnames(mat), sppOut))
    })
    
    ## return scale as area, and metric summary
    out <- do.call(rbind, out)
    
    return(out)
}

## helper function to determine scale for raster of plot
## note: scale (or resolution) must be a multiple of each dimension of the plot
## to work easily, so just scale each side independently
#' @param r is the raster object of the plot produced internally in `scaleMetric`

determineScale <- function(r) {
    xmx <- xmax(r)
    ymx <- ymax(r)
    xIsMax <- xmx >= ymx
    
    ## target minimum area
    targetArea <- 10
    
    ## number of halvings allowed by dimensions of plot
    n <- floor((log(xmx * ymx) - log(targetArea)) / log(2) / 2)
    
    s <- rep(0:n, each = 2)
    
    if(xIsMax) {
        return(cbind(xmx / 2^c(0, rep(1:n, each = 2)),
                     ymx / 2^rep(0:n, each = 2)[-(n+1)*2]))
    } else {
        return(cbind(xmx / 2^rep(0:n, each = 2)[-(n+1)*2],
                     ymx / 2^c(0, rep(1:n, each = 2))))
    }
}


## calculate singletons across scales
negbBCI <- scaleNegBinom(bci)
negbPASO <- scaleNegBinom(paso)
