## =======================================
## functions to aid in scaling SAD metrics
## =======================================

library(pika)
library(raster)

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
    # browser()
    
    ## loop over scales, cutting data by re-scaled raster and calculating SAD z-values
    out <- lapply(1:nrow(scales), function(s) {
        ## rescale
        res(r) <- scales[s, ]
        print(s)
        
        ## within a scale, loop over cells to calcuate SAD z-values within each
        
        cells <- cellFromXY(r, xy = x[, c('x', 'y')])
        cellsPerm <- sample(cells) # permuted cells
        
        ## make sure we're never sampling more than 2^(10 - 4) == 64 cells
        if(s > 4) {
            theseCells <- sample(unique(cells), 64)
        } else {
            theseCells <- unique(cells)
        }
        
        allM <- lapply(theseCells, FUN = function(i) {
            m <- calcAtCell(x, cells, i, fun)
            
            if(s == 1 | !perm) {
                mPerm <- m
            } else {
                mPerm <- calcAtCell(x, cellsPerm, i, fun)
            }
            
            return(c(m, perm = mPerm))
        })
        
        allM <- do.call(rbind, allM)
        
        ## return scale and metric values
        return(cbind(scale = prod(scales[s, ]), allM))
    })
    
    ## return scale as area, and metric summary
    out <- as.data.frame(do.call(rbind, out))
    
    return(out)
}


## helper function to subset data by cell and calculate metric within that cell
#' @param x is the data
#' @param cells is a matrix of cell IDs
#' @param i is the row in `cells` indicating the cell of interest
#' @param fun is the function to calculate the metric

calcAtCell <- function(x, cells, i, fun) {
    newx <- x[cells == i, ]
    
    return(fun(x, newx))
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
