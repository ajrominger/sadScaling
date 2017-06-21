## =======================================
## functions to aid in scaling SAD metrics
## =======================================

library(pika)
library(raster)

## function to take plot data, divide it up across scales, and caluclate
## metrics across scales
#' @param x is the data
#' @param fun is a function returning the metric of interest
#' @param mod is the model to be fit
#' @param perm is a logical indicating whether to compute the metric on permuted cells

scaleMetric <- function(x, fun, mod = 'fish', perm = TRUE) {
    ## make raster object for plot
    r <- raster(ncols = ceiling(max(x$x)), nrows = ceiling(max(x$y)), 
                xmn = 0, xmx = ceiling(max(x$x)), 
                ymn = 0, ymx = ceiling(max(x$y)))
    
    ## determine range of scales
    scales <- determineScale(r)
    
    ## loop over scales, cutting data by re-scaled raster and calculating SAD z-values
    out <- mclapply(1:nrow(scales), mc.cores = 6, FUN = function(s) {
        ## rescale
        res(r) <- scales[s, ]
        
        ## within a scale, loop over cells to calcuate SAD z-values within each
        
        cells <- cellFromXY(r, xy = x[, c('x', 'y')])
        cellsPerm <- sample(cells) # permuted cells
        
        allM <- lapply(unique(cells), FUN = function(i) {
            m <- calcAtCell(x, cells, i)
            
            if(s == 1 | !perm) {
                mPerm <- m
            } else {
                mPerm <- calcAtCell(x, cellsPerm, i)
            }
            
            return(c(m, mPerm))
        })
        
        allM <- do.call(rbind, ll)
        
        ## return averaged metric vales
        return(c(mean = mean(allM[, 1]), 
                 ci = quantile(allM[, 1], c(0, 0.95), names = FALSE),
                 mPermMean = mean(allM[, 2]), 
                 mPermCI = quantile(allM[, 2], c(0, 0.95), names = FALSE)))
    })
    
    ## return scale as area, and metric summary
    out <- data.frame(scale = scales[, 1] * scales[, 2], do.call(rbind, out))
    
    return(out)
}


## helper function to subset data by cell and calculate metric within that cell
#' @param x is the data
#' @param cells is a matrix of cell IDs
#' @param i is the row in `cells` indicating the cell of interest

calcAtCell <- function(x, cells, i) {
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
#' @param r is the raster object of the plot produced internally in `scaleMetric`

determineScale <- function(r) {
    xmx <- xmax(r)
    ymx <- ymax(r)
    
    ## target minimum area
    targetArea <- 10
    
    ## number of halvings allowed by dimensions of plot
    n <- floor(log(targetArea / (xmx * ymx)) / (-2 * log(2)))
    
    return(cbind(xmx * 2^-(0:n), ymx * 2^-(0:n)))
}
