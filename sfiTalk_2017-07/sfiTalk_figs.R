## script for talk at SFI postdoc conference July 2017

library(igraph)
library(ape)
library(socorro)
library(pika)

setwd('~/Dropbox/Research/sadScaling/sfiTalk_2017-07')

## function implementing the foodweb niche model
nicheMod <- function(S,K,ntry=1) {
    if(ntry < 10) {
        # y <- 1 + (1-2*K)*log(2)/log(K)
        # x <- 2*K
        y <- 1
        x <- 1
        
        bet <- 1/(2*K^y) - 1
        
        ni <- runif(S)
        
        ri <- rbeta(S,1,bet)*ni^x
        ri[which.min(ni)] <- 0
        
        ci <- runif(S,ri/2,ni)
        
        feed.min <- ci - ri/2
        feed.max <- ci + ri/2
        
        dirct <- mapply(FUN = function(x1,x2) {
            1*(ni <= x2 & ni >= x1)
        },
        feed.min,feed.max)
        
        undirct <- dirct + t(dirct)
        undirct[undirct > 1] <- 1
        
        all0 <- which(rowSums(undirct)==0)
        undirct <- undirct[-all0,-all0]
        
        if(nrow(undirct) == 0) {
            nicheMod(S,K,ntry+1)
        } else {
            return(undirct)
        }
    } else {
        out <- matrix(0,nrow=S-1,ncol=S-1)
        diag(out) <- 1
        out <- cbind(rbind(0,out),0)
        out <- out + t(out)
        out[out > 1] <- 1
        return(out)
    }
}



set.seed(1)
x <- nicheMod(100, 0.04)
x <- x[rowSums(x) > 0, colSums(x) > 0]
x[lower.tri(x, diag = TRUE)] <- 0
N <- nrow(x)
x <- graph_from_adjacency_matrix(x, mode = 'undirected')

pdf('fig_trophicWeb.pdf', width = 4, height = 4)
lay <- layout_as_tree(x)
lay[lay[, 2] == min(lay[, 2]), 2] <- min(lay[lay[, 2] != min(lay[, 2]), 2])
par(mar = rep(0.1, 4), bg = 'black')
plot(x, vertex.label = NA, vertex.size = 2, 
     vertex.color = 'white', vertex.frame.color = 'white',
     edge.color = hsv(0, 0, 1, alpha = 0.25),
     layout = lay)
dev.off()

set.seed(2)
tre <- rphylo(N, 0.9, 0.5, fossils = TRUE)
e <- tre$edge[, 2]
eExtant <- which(e[e %in% 1:Ntip(tre)] %in% 1:N)
lay2 <- cbind(eExtant, rep(1, N))

pdf('fig_trophicWeb-phylo.pdf', width = 6, height = 4)
layout(matrix(2:1, nrow = 2), heights = c(1, 3))
par(bg = 'black', fg = 'white', mar = rep(0, 4))

plot(tre, direction = 'upward', edge.color = 'white', show.tip.label = FALSE, 
     y.lim = c(-0.3, 7.5), yaxs = 'i')

plot(x, vertex.label = NA, vertex.size = 200, 
     vertex.color = 'white', vertex.frame.color = 'white',
     edge.color = hsv(0, 0, 1, alpha = 0.25), edge.curved = 1,
     layout = lay2, ylim = c(1, 70), xlim = c(1, Ntip(tre)), 
     rescale = FALSE, asp = 0)

dev.off()


set.seed(1)
x <- nicheMod(100, 0.4)
x <- x[rowSums(x) > 0, colSums(x) > 0]
N <- nrow(x)
x <- x * matrix(runif(N^2), nrow = N)

image(x, asp = 1, col = colorRampPalette(hsv(c(0.17, 0.20, 0.30, 0.60, 0.80), 
                                             c(0.30, 0.80, 1.00, 0.80, 0.80), 
                                             c(0.90, 0.80, 1.00, 0.80, 0.80)))(20)))

pdf('fig_trophicMat.pdf', width = 4, height = 4)
lay <- layout_as_tree(x)
lay[lay[, 2] == min(lay[, 2]), 2] <- min(lay[lay[, 2] != min(lay[, 2]), 2])
par(mar = rep(0.1, 4), bg = 'black')
plot(x, vertex.label = NA, vertex.size = 2, 
     vertex.color = 'white', vertex.frame.color = 'white',
     edge.color = hsv(0, 0, 1, alpha = 0.25),
     layout = lay)
dev.off()

set.seed(2)
tre <- rphylo(N, 0.9, 0.5, fossils = TRUE)
e <- tre$edge[, 2]
eExtant <- which(e[e %in% 1:Ntip(tre)] %in% 1:N)
lay2 <- cbind(eExtant, rep(1, N))

pdf('fig_trophicWeb-phylo.pdf', width = 6, height = 4)
layout(matrix(2:1, nrow = 2), heights = c(1, 3))
par(bg = 'black', fg = 'white', mar = rep(0, 4))

plot(tre, direction = 'upward', edge.color = 'white', show.tip.label = FALSE, 
     y.lim = c(-0.3, 7.5), yaxs = 'i')

plot(x, vertex.label = NA, vertex.size = 200, 
     vertex.color = 'white', vertex.frame.color = 'white',
     edge.color = hsv(0, 0, 1, alpha = 0.25), edge.curved = 1,
     layout = lay2, ylim = c(1, 70), xlim = c(1, Ntip(tre)), 
     rescale = FALSE, asp = 0)

dev.off()

## environmental variables

set.seed(3)
f1 <- spline(runif(10), cumsum(rnorm(10)), n = 100)
f2 <- spline(runif(10), cumsum(rnorm(10)), n = 100)

pdf('fig_trophicWeb-phylo.pdf', width = 1.5, height = 3)
par(bg = 'black', fg = 'white', col.lab = 'white', mfrow = 1:2, 
    mar = c(2, 0, 0, 0) + 0.2, mgp = c(1.75, 0.5, 0))

plot(f1$y, f1$x, 
     type = 'l', lwd = 3, col = hsv(0.2, 0.5, 1), 
     xlab = '', ylab = '', axes = FALSE, frame.plot = TRUE)
mtext('Env1', side = 1, line = 0.5)

plot(f2$y, f2$x, 
     type = 'l', lwd = 3, col = hsv(0.1, 0.5, 1), 
     xlab = '', ylab = '', axes = FALSE, frame.plot = TRUE)
mtext('Env2', side = 1, line = 0.5)

dev.off()


## sad motivation and theory fig