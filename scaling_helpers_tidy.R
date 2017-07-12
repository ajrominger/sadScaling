## =======================================
## functions to aid in scaling SAD metrics
## =======================================

library(pika)
library(raster)
library(parallel)

## function to take plot data, divide it up across scales, and caluclate
## metrics across scales
#' @param x is the data
#' @param fun is a function returning the metric of interest
#' @param perm is a logical indicating whether to compute the metric on permuted cells

scaleMetric <- function(x, fun, perm = TRUE) {
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
        
        ## make sure we're never sampling more than 2^5 cells
        if(s > 6) {
            theseCells <- sample(unique(cells), 2^5)
        } else {
            theseCells <- unique(cells)
        }
        
        allM <- lapply(theseCells, FUN = function(i) {
            newx <- x[cells == i, ]
            return(c(N = sum(newx$count), fun(x, newx)))
        })
        
        allM <- data.frame(type = 'spatial', do.call(rbind, allM))
        Nbar <- round(mean(allM$N))
        
        if(s != 1) {
            permM <- lapply(1:2^5, function(i) {
                newx <- x[sample(nrow(x), Nbar), ]
                return(c(N = Nbar, fun(x, newx)))
            })
            permM <- data.frame(type = 'random', do.call(rbind, permM))
        } else {
            permM <- allM
        }
        
        ## return scale and metric values
        return(cbind(scale = prod(scales[s, ]), rbind(allM, permM)))
    })
    
    ## return scale as area, and metric
    out <- as.data.frame(do.call(rbind, out))
    
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
