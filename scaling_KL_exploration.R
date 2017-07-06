## script exploring scaling of KL distance (between NB and pois) when the mean 
## of the distrib increases.  shows that KL gets bigger with increasing mean,
## even if other params of the distribution generating the data are the same.
## so scaling of KL distance shouldn't be considered too seriously because 
## mean will increase with area

nout <- 6
M <- seq(1, 100, length.out = nout)


layout(matrix(c(rep(nout + 1, nout), 1:nout), nrow = 2, byrow = TRUE))
par(mar = rep(0.1, 4), oma = c(2, 2, 0, 0))
kl <- sapply(M, function(mu) {
    x <- rnbinom(1000, size = 200, mu = mu)
    nbpar <- fitdistr(x, 'negative binomial')$estimate
    poispar <- mean(x)
    p <- dnbinom(0:10^6, size = nbpar[1], mu = nbpar[2], log = TRUE)
    q <- dpois(0:10^6, poispar, log = TRUE)
    
    plot(0:120, exp(p[1:121]), type = 'l', col = 'blue', ylim = exp(range(p[1:121], q[1:121])))
    points(0:120, exp(q[1:121]), type = 'l', col = 'red')
    
    sum(exp(p) * (p - q))
})

plot(M, kl, type = 'b')
