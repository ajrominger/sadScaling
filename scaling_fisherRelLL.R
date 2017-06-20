library(pika)
library(socorro)
library(raster)
library(parallel)

setwd('~/Dropbox/Research/sadScaling')

source('../happySAD/sampleTransform.R')

bci <- read.csv('../data/stri/BCIS.csv', as.is = TRUE)
bci <- bci[bci$year == max(bci$year), ]


## function to take plot data, divide it up across scales, and caluclate
## summaries of relative log likelihoods across scale
scaleRelLL <- function(x) {
    ## make raster object for plot
    r <- raster(ncols = ceiling(max(x$x)), nrows = ceiling(max(x$y)), 
                xmn = 0, xmx = ceiling(max(x$x)), 
                ymn = 0, ymx = ceiling(max(x$y)))
    
    ## determine range of scales
    scales <- determineScale(r)
    
    ## loop over scales, cutting data by re-scaled raster and calculating SAD z-values
    out <- mclapply(1:nrow(scales), mc.cores = 6, FUN = function(s) {
        if(s == 1) {
            abund <- tapply(x$count, x$spp, sum)
            thisSAD <- sad(abund, 'fish')
            st <- sampTrans(abund, sad(abund), sum(abund), include0 = FALSE, log = TRUE)
            
            return(c(llMean = logLik(thisSAD) - sum(st), llCI = rep(NA, 2),
                     llPermMean = NA, llPermCI = rep(NA, 2)))
        } else {
            ## rescale
            res(r) <- scales[s, ]
            
            ## within a scale, loop over cells to calcuate SAD z-values within each
            
            cells <- cellFromXY(r, xy = x[, c('x', 'y')])
            # cells[is.na(cells)] <- -1 # avoid issues with NA when matching below
            cellsPerm <- sample(cells) # permuted cells
            
            ll <- lapply(unique(cells), FUN = function(i) {
                ll <- calcLL(x, cells, i)
                llPerm <- calcLL(x, cellsPerm, i)
                return(c(ll, llPerm))
            })
            
            ll <- do.call(rbind, ll)
            
            ## return averaged z vales
            return(c(llMean = mean(ll[, 1] - ll[, 2]), 
                     llCI = quantile(ll[, 1] - ll[, 2], c(0, 0.95), names = FALSE),
                     llPermMean = mean(ll[, 3] - ll[, 4]), 
                     llPermCI = quantile(ll[, 3] - ll[, 4], c(0, 0.95), names = FALSE)))
        }
    })
    
    ## return scale as area, and z-value summary
    out <- data.frame(scale = scales[, 1] * scales[, 2], do.call(rbind, out))
    
    return(out)
}

## helper function to subset data by cell and calculate log liks
calcLL <- function(x, cells, i) {
    newx <- x[cells == i, ]
    abund <- tapply(newx$count, newx$spp, sum)
    llF <- logLik(sad(abund, 'fish'))
    llST <- sum(sampTrans(abund, sad(tapply(x$count, x$spp, sum)), nrow(newx), 
                          include0 = FALSE, log = TRUE))
    
    return(c(llF, llST))
}


## helper function to determine scale for raster of plot
## note: scale (or resolution) must be a multiple of each dimension of the plot
## to work easily, so just scale each side independently
determineScale <- function(r) {
    xmx <- xmax(r)
    ymx <- ymax(r)
    
    ## target minimum area
    targetArea <- 10
    
    ## number of halvings allowed by dimensions of plot
    n <- floor(log(targetArea / (xmx * ymx)) / (-2 * log(2)))
    
    return(cbind(xmx * 2^-(0:n), ymx * 2^-(0:n)))
}

bciScaleRelLL <- scaleRelLL(bci)
write.csv(bciScaleRelLL, 'scaling_fisherRelLL_BCI.csv')
