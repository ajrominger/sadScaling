library(pika)
library(parallel)
library(socorro)
library(sp)
library(raster)

setwd('~/Dropbox/Research/sadScaling')
dataWD <- '~/Dropbox/Research/data/stri'

logLikZ_old <- function(x, nrep=1000, return.sim=FALSE) {
    lik.obs <- logLik(x)
    n <- x$nobs
    rfun <- getrfun(x)
    dfun <- getdfun(x)
    
    lik.sim <- replicate(nrep, {
        newx <- rfun(n)
        sum(dfun(newx, log=TRUE))
    })
    
    z <- ((lik.obs - mean(lik.sim)) / sd(lik.sim))^2
    
    if(return.sim) {
        lik.sim <- ((lik.sim - mean(lik.sim)) / sd(lik.sim))^2
    } else {
        lik.sim <- NULL
    }
    
    return(list(z=as.numeric(z), obs=lik.obs, sim=lik.sim))
}

# x <- read.csv(file.path(dataWD, 'PASO.csv'), as.is = TRUE)
# x <- read.csv(file.path(dataWD, 'BCIS.csv'), as.is = TRUE)
x <- read.csv(file.path(dataWD, 'KORU.csv'), as.is = TRUE)
x <- x[x$year == max(x$year), ]

N <- exp(seq(log(100), log(nrow(x)), length.out = 10))

sadSumm <- lapply(N, function(n) {
    out <- replicate(10, {
        newx <- x[sample(1:nrow(x), n), ]
        abund <- tapply(newx$count, newx$spp, sum)
        n1 <- sum(abund == 1)
        s <- length(abund)
        z <- logLikZ(sad(abund, 'fish', keepData = TRUE))$z
        
        return(c(n1 = n1, s = s, z = z))
    })
    
    return(rowMeans(out))
})

sadSumm <- data.frame(N = N, do.call(rbind, sadSumm))

plot(sadSumm$N, sadSumm$z, log = 'xy')
plot(sadSumm$s, sadSumm$z)
plot(sadSumm$N, sadSumm$n1, log = 'xy')

plot(sadSumm$n1, sadSumm$s, col = quantCol(sadSumm$z, 
                                           pal = c('yellow', 'cyan', 'blue'), 
                                           trans = 'log'))
plot(sadSumm$N, sadSumm$s, col = quantCol(sadSumm$z, 
                                          pal = c('yellow', 'cyan', 'blue'), 
                                          trans = 'log'))

plot(sadSumm$N, sadSumm$z, log = 'xy')
