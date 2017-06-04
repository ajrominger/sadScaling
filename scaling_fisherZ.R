library(pika)
library(parallel)
library(socorro)
library(sp)
library(raster)

setwd('~/Dropbox/Research/sadScaling')

## source hidden function from meteR to deal with areas
source('~/Dropbox/Research/meteR/R/sar_helper_funs.R')

## wd storing the data
dataWD <- '~/Dropbox/Research/data/stri'


## function to take plot data, divide it up across scales, and caluclate
## summaries of z-values across scale
scaleZ <- function(x) {
    ## make raster object for plot
    r <- raster(ncols = ceiling(max(x$x)), nrows = ceiling(max(x$y)), 
                xmn = 0, xmx = ceiling(max(x$x)), 
                ymn = 0, ymx = ceiling(max(x$y)))
    
    ## determine range of scales
    scales <- determineScale(r)
    
    ## loop over scales, cutting data by re-scaled raster and calculating SAD z-values
    out <- lapply(1:nrow(scales), function(s) {
        if(s == 1) {
            abund <- tapply(x$count, x$spp, sum)
            thisSAD <- sad(abund, 'fish', keepData = TRUE)

            return(c(zMean = logLikZ(thisSAD)$z, zCI = rep(NA, 2),
                     zPermMean = NA, zPermCI = rep(NA, 2)))
        } else {
            ## rescale
            res(r) <- scales[s, ]

            ## within a scale, loop over cells to calcuate SAD z-values within each

            cells <- cellFromXY(r, xy = x[, c('x', 'y')])
            # cells[is.na(cells)] <- -1 # avoid issues with NA when matching below
            cellsPerm <- sample(cells) # permuted cells

            # zz <- mclapply(unique(cells[cells > 0]), mc.cores = 6, function(i) {
            zz <- mclapply(unique(cells), mc.cores = 6, function(i) {
                newx <- x[cells == i, ]
                abund <- tapply(newx$count, newx$spp, sum)
                thisSAD <- sad(abund, 'fish', keepData = TRUE)
                z <- calcZ(x, cells, i)
                zPerm <- calcZ(x, cellsPerm, i)
                return(c(z, zPerm))
            })

            zz <- do.call(rbind, zz)


            ## return averaged z vales
            return(c(zMean = mean(zz[, 1]), zCI = quantile(zz[, 1], c(0, 0.95),
                                                           names = FALSE),
                     zPermMean = mean(zz[, 2]), zPermCI = quantile(zz[, 2], c(0, 0.95),
                                                                   names = FALSE)))
        }
    })

    ## return scale as area, and z-value summary
    out <- data.frame(scale = scales[, 1] * scales[, 2], do.call(rbind, out))

    return(out)
}

## helper function to subset data by cell and calculate z-vals
calcZ <- function(x, cells, i) {
    newx <- x[cells == i, ]
    abund <- tapply(newx$count, newx$spp, sum)
    thisSAD <- sad(abund, 'fish', keepData = TRUE)
    return(logLikZ(thisSAD)$z)
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


## loop over data, computing z scaling for all
allZ <- lapply(list.files(dataWD), function(f) {
    ## read in data
    x <- read.csv(file.path(dataWD, f), as.is = TRUE)
    
    ## get most recent census
    x <- x[x$year == max(x$year), ]
    
    ## get data.frame of z value scaling
    out <- cbind(plot = gsub('.csv', '', f), scaleZ(x))
    return(out)
})

allZ <- do.call(rbind, allZ)

write.csv(allZ, 'scaling_fisherZ.csv', row.names = FALSE)
