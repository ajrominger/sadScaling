library(pika)
library(parallel)
library(socorro)
library(sp)
library(raster)

setwd('~/Dropbox/Research/maxentConcept')

# ## source hidden function from meteR to deal with areas
# source('~/Dropbox/Research/meteR/R/sar_helper_funs.R')

## wd storing the data
dataWD <- '~/Dropbox/Research/data/stri'

## read in data
x <- read.csv(file.path(dataWD, list.files(dataWD)[1]), as.is = TRUE)

## get most recent census
x <- x[x$year == max(x$year), ]

## function to take plot data, divide it up across scales, and caluclate
## summaries of z-values across scale

scaleZ <- function(x) {
    ## make raster object for plot
    r <- raster(ncols = ceiling(max(x$x)), nrows = ceiling(max(x$y)), 
                xmn = 0, xmx = ceiling(max(x$x)), 
                ymn = 0, ymx = ceiling(max(x$y)))
    
    ## determine range of scales
    scales <- 2^(2:floor(log(min(max(xmax(r), ymax(r)) / 2, 
                                 min(xmax(r), ymax(r))), base = 2)))
    
    ## loop over scales, cutting data by re-scaled raster and calculating SAD z-values
    out <- lapply(scales, function(s) {
        res(r) <- s
        
        ## within a scale, loop over cells to calcuate SAD z-values within each
        
        cells <- cellFromXY(r, xy = x[, c('x', 'y')])
        cells[is.na(cells)] <- -1 # avoid issues with NA when matching below
        cellsPerm <- sample(cells) # permuted cells
        
        zz <- mclapply(unique(cells[cells > 0]), mc.cores = 6, function(i) {
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
    })
    
    out <- cbind(scale = scales, do.call(rbind, out))
    
    return(out)
}

## helper function to subset data by cell and calculate z-vals
calcZ <- function(x, cells, i) {
    newx <- x[cells == i, ]
    abund <- tapply(newx$count, newx$spp, sum)
    thisSAD <- sad(abund, 'fish', keepData = TRUE)
    return(logLikZ(thisSAD)$z)
}

foo <- scaleZ(x)
