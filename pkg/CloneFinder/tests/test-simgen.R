library(CloneFinder)

pars.default <- list('rounds'=400, 'nu'=0, 'pcnv'=1, 'norm.contam'=FALSE)
dataPars <- list('snps.seq'=1000000,'snps.cgh'=600000,'mu'=70,'sigma.reads'=25,'sigma0.lrr'=.15, 
                 'sigma0.baf'=.03,'density.sigma'=.1)
threshold <- .04

psis <- c(0.6, 0.3, 0.1)

tumor <- tumorGen(psis, pars.default$rounds, pars.default$nu, pars.default$pcnv, pars.default$norm.contam)

data <- dataGen(tumor, dataPars$snps.seq, dataPars$snps.cgh, dataPars$mu, dataPars$sigma.reads, 
                dataPars$sigma0.lrr, dataPars$sigma0.baf, dataPars$density.sigma)



plot.data(tumor, data)
