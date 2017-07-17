library(socorro)
library(plyr)

setwd('~/Dropbox/Research/sadScaling')

negbScale <- read.csv('scaling_negbiomSSAD.csv', as.is = TRUE)

negbScaleLLSumm <- ddply(negbScale, c('site', 'scale'), function(x) {
    d <- x$nbLL - x$pLL
    permD <- x$perm.nbLL - x$perm.pLL
    
    return(c(meanD = mean(d), ciD = quantile(d, c(0.025, 0.975), names = FALSE), 
             meanPermD = mean(permD), 
             ciPermD = quantile(permD, c(0.025, 0.975), names = FALSE)))
})



bci <- read.csv('../data/stri/BCIS.csv', as.is = TRUE)
bci <- bci[bci$year == max(bci$year), ]

paso <- read.csv('../data/stri/PASO.csv', as.is = TRUE)
paso <- paso[paso$year == max(paso$year), ]





foo <- tidy2mat(negbScale$scale[negbScale$site == 'BCI'], negbScale$spp[negbScale$site == 'BCI'],
                negbScale$nbPar.k[negbScale$site == 'BCI'])

foofa <- psych::factor.pa(foo, nfactors = 6)


bciNB <- negbScale[negbScale$site == 'BCI', ]
bciNB$f <- apply(unclass(foofa$loadings), 1, which.max)


bla <- ddply(bciNB, c('scale', 'f'), function(x) c(k = mean(x$nbPar.k), perm.k = mean(x$perm.nbPar.k)))

plot(1, xlim = range(bla$scale), ylim = range(bla[, 3:4]), log = 'xy', type = 'n')
for(i in sort(unique(bla$f))) {
    points(bla$scale[bla$f == i], bla$k[bla$f == i], col = hsv(i/10), type = 'l')
}
abline(h = 100)




plot(1, xlim = range(bciNB$scale), ylim = range(bciNB$nbPar.k), type = 'n', log = 'xy')
for(i in sort(unique(f))) {
    d_ply(negbScale[negbScale$site == 'BCI' & negbScale$spp %in% names(which(f == i)), ], 'spp', function(x) {
        points(x$scale, x$nbPar.k, type = 'l', col = hsv(i / 10, alpha = 0.5))
    })
}

plot(1, xlim = range(negbScale$scale), ylim = range(negbScale$nbPar.k), type = 'n', log = 'xy')
d_ply(negbScale[negbScale$site == 'BCI', ], 'spp', function(x) {
    points(x$scale, x$nbPar.k, type = 'l', col = gray(0.5, alpha = 0.5))
})
