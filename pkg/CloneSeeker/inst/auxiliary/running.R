library(CloneFinder)
library(foreach)
library(doParallel)
#KRC: source('parOpt/clonefinder/functions.R')
#KRC: source('parOpt/clonefinder/algorithm.R')

convert.sc <- function(mutdata, cndata){
  segs <- mutdata$seg
  segs.unique <- unique(segs)
  sciclone.cn <- as.data.frame(t(sapply(1:nrow(cndata), function(k){
    start <- sum(cndata$markers[1:k])+1 - cndata$markers[k]
    end <- sum(cndata$markers[1:k])
    c('chr'=cndata$chr[k], 'start'=start, 'stop'=end, 'segment_mean'=cndata$X[k] + cndata$Y[k])
  })))
  #randomly generate loci for mutations on their respective segments:
  if(length(segs.unique)>0){
    mutcn <- round(t(sapply(1:nrow(mutdata), function(k){
      sort(c(cndata$X[segs[k]], cndata$Y[segs[k]]), decreasing=TRUE)
    })))
    colnames(mutcn) <- c('A', 'B')
    genotypes <- sapply(1:nrow(mutcn),function(k){
      paste(c(rep(colnames(mutcn)[1], mutcn[k,1]), rep(colnames(mutcn)[2], mutcn[k,2])), collapse='')
    })
    loci <- unlist(lapply(1:length(segs.unique),function(k){
      minim <- sciclone.cn$start[segs.unique[k]]
      maxim <- sciclone.cn$stop[segs.unique[k]]
      n <- length(which(segs==segs.unique[k]))
      sample(minim:maxim, n, replace=FALSE)
    }))
  }else{
    loci <- numeric(0)
  }
  sciclone.dat <- data.frame('chr'=mutdata$chr, 'start'=loci, 'refCounts'=mutdata$refCounts, 
                             'varCounts'=mutdata$varCounts, 'VAF'=mutdata$VAF, 'mutid'=mutdata$mut.id)
  list('sciclone.dat'=sciclone.dat, 'sciclone.cn'=sciclone.cn)
}

convert.exp <- function(mutdata, cndata){
  segs <- mutdata$seg
  segs.unique <- unique(segs)
  start <- sapply(1:nrow(cndata),function(k){sum(cndata$markers[1:k])+1 - cndata$markers[k]})
  end <- sapply(1:nrow(cndata),function(k){sum(cndata$markers[1:k])})
  CBS <- data.frame('chr'=cndata$chr, 'startpos'=start, 'endpos'=end, 'CN_Estimate'=cndata$X+cndata$Y)
  if(length(segs.unique)>0){
    loci <- unlist(lapply(1:length(segs.unique),function(k){
      minim <- CBS$startpos[segs.unique[k]]
      maxim <- CBS$endpos[segs.unique[k]]
      n <- length(which(segs==segs.unique[k]))
      sample(minim:maxim, n, replace=FALSE)
    }))
  }else{
    loci <- numeric(0)
  }
  mutids <- mutdata$mut.id
  SNV <- data.frame('chr'=mutdata$chr, 'startpos'=loci, 'AF_Tumor'=mutdata$VAF, 'PN_B'=rep(0, nrow(mutdata)), 'mutid'=mutids)
  list('SNV'=as.matrix(SNV), 'CBS'=as.matrix(CBS))
}

runAlgs <- function(i, set, minim.vaf, datapath, respath, mu, cf=TRUE, sc=FALSE, exp=TRUE, cn=TRUE, mix, real.dat, run, pars){
  datafiles <- list.files(datapath)
  dat <- get(load(paste(datapath, '/', datafiles[i], sep='')))
  number <- strsplit(strsplit(datafiles[i], split='-')[[1]][3], split='.rda')[[1]][1]
  if(mix | real.dat){
    number <- strsplit(strsplit(datafiles[i], split='-')[[1]][2], split='.rda')[[1]][1]
  }
  if(length(dat)==1){
    vardata <- dat$dat$seq.data
    cndata <- dat$dat$cn.data
  }else{
    vardata <- dat$seq.data
    cndata <- dat$cn.data
  }
  mutdata <- vardata[vardata$status=='somatic',]
  mut.filt <- filter.mut(mutdata, mu=mu, threshold=3)
  mutdata.filt <- mut.filt$mat
  mutdata.filt <- mutdata.filt[mutdata.filt$chr %in% 1:22,]
  mutdata.filt <- mutdata.filt[which(mutdata.filt$VAF > minim.vaf/2),]
  if(cn){
    if(cf){
      t1 <- Sys.time()
      res.cf <- runAlg(cndata, vardata, cnmodels, psisets[[2]], pars=pars, priorbank=NULL)
      t2 <- Sys.time()
      runtime.cf <- difftime(t2, t1, units='secs')
      res.cf <- list('res'=res.cf, 'runtime'=runtime.cf)
      if(!is.null(run)){
        filename <- paste(respath, '/res-cf-', set, run, '/res-cf-', set, '-',  number, '.rda', sep='')
      }else{
        filename <- paste(respath, '/res-cf-', set, '/res-cf-', set, '-',  number, '.rda', sep='')
      }
      save(res.cf, file=filename)
    }
    if(exp){
      inputs.exp <- convert.exp(mutdata.filt, cndata)
      t1 <- Sys.time()
      res.exp <- runExPANdS(inputs.exp$SNV,inputs.exp$CBS, precision=.04, min_CellFreq = minim.vaf)
      t2 <- Sys.time()
      runtime.exp <- difftime(t2, t1, units='secs')
      res.exp <- list('res'=res.exp, 'runtime'=runtime.exp)
      filename <- paste(respath, '\\res-exp-', set,'\\res-exp-', set, '-',  number, '.rda', sep='')
      save(res.exp, file=filename)
    }
    if(sc){
      inputs.sc <- convert.sc(mutdata.filt, cndata)
      t1 <- Sys.time()
      res.sc <- sciClone(inputs.sc$sciclone.dat, inputs.sc$sciclone.cn, sampleNames=paste('sample-',i, sep=''), maximumClusters = 10, minimumDepth = 50)
      t2 <- Sys.time()
      runtime.sc <- difftime(t2, t1, units='secs')
      res.sc <- list('res'=res.sc, 'runtime'=runtime.sc)
      filename <- paste(respath, '/res-sc-', set,'/res-sc-', set, '-', number, '.rda', sep='')
      save(res.sc, file=filename)
    }
  }else{
    if(cf){
      t1 <- Sys.time()
      res.cf <- runAlg(cndata, vardata, cnmodels, psisets[[2]], pars=pars, priorbanks=NULL)
      t2 <- Sys.time()
      runtime.cf <- difftime(t2, t1, units='secs')
      res.cf <- list('res'=res.cf, 'runtime'=runtime.cf, 'parameters'=c(list('mu'=mu,'minim.vaf'=minim.vaf),pars))
      filename <- paste(respath, '/res-cf-', set, '/res-cf-', set, '-',  number, '.rda', sep='')
      save(res.cf, file=filename)
    }
    if(exp){
      inputs.exp <- convert.exp(mutdata.filt, cndata)
      t1 <- Sys.time()
      res.exp <- res.exp <- runExPANdS(inputs.exp$SNV, precision=.04, min_CellFreq = .05)
      t2 <- Sys.time()
      runtime.exp <- difftime(t2, t1, units='secs')
      res.exp <- list('res'=res.exp, 'runtime'=runtime.exp)
      filename <- paste(respath, '\\res-exp-', set,'\\res-exp-', set, '-',  number, '.rda', sep='')
      save(res.exp, file=filename)
    }
    if(sc){
      inputs.sc <- convert.sc(mutdata.filt, cndata)
      t1 <- Sys.time()
      res.sc <- sciClone(inputs.sc$sciclone.dat, sampleNames=paste('sample-',i, sep=''), maximumClusters = 10, minimumDepth = 50)
      t2 <- Sys.time()
      runtime.sc <- difftime(t2, t1, units='secs')
      res.sc <- list('res'=res.sc, 'runtime'=runtime.sc)
      filename <- paste(respath, '/res-sc-', set,'/res-sc-', set, '-', number, '.rda', sep='')
      save(res.sc, file=filename)
    }
  }
}


###SciClone and Clonefinder, on the cluster:
#Loading and defining stuff:
kmax <- 5
load('parOpt/clonefinder/cnmodels.rda')
load('parOpt/clonefinder/psis.10.rda')
load('parOpt/clonefinder/psis.20.rda')
load('parOpt/clonefinder/psis.50.rda')
load('parOpt/clonefinder/psis.100.rda')
load('parOpt/clonefinder/etamat.10.rda')
load('parOpt/clonefinder/etamat.20.rda')
load('parOpt/clonefinder/etamat.50.rda')
psisets <- list(psis.10, psis.20,psis.50,psis.100)
etasets <- list(etamat.10, etamat.20,etamat.50)

setNumber <- 'CLL'
dpath <- paste('parOpt/dat/data-', setNumber, sep='')
rpath <- paste('parOpt/res/res-', setNumber, sep='')
#metadata <- get(load(paste('parOpt/sim/sims-',setNumber,'/metadata-',setNumber,'.rda', sep='')))
cf.var <- TRUE
#mix <- real.dat <- sc.var <- FALSE
cf.var <- real.dat <- TRUE
mix <- sc.var <- FALSE
run <- 5
X <- list(1:262)
nodes <- 28
parlist <- list(list('sigma0'=.4, 'theta'=.5, 'ktheta'=.3, 'mtheta'=.5, 'alpha'=.5, 'thresh'=0, 'cutoff'=40, 'Q'=100, 'iters'=4),
                list('sigma0'=.3, 'theta'=.5, 'ktheta'=.3, 'mtheta'=.5, 'alpha'=.5, 'thresh'=0, 'cutoff'=40, 'Q'=100, 'iters'=4),
                list('sigma0'=8, 'theta'=.5, 'ktheta'=.5, 'mtheta'=.5, 'alpha'=.5, 'thresh'=.04, 'cutoff'=500, 'Q'=100, 'iters'=4),
                list('sigma0'=6, 'theta'=.6, 'ktheta'=.6, 'mtheta'=.5, 'alpha'=.5, 'thresh'=.04, 'cutoff'=500, 'Q'=100, 'iters'=4),
                list('sigma0'=5, 'theta'=.9, 'ktheta'=.3, 'mtheta'=.9, 'alpha'=.5, 'thresh'=.04, 'cutoff'=100, 'Q'=100, 'iters'=4),
                list('sigma.lrr'=.15,'sigma.baf'=.03, 'sigmaM'=.14, 'theta'=.9, 'ktheta'=.9, 'mtheta'=.9, 'alpha'=.5, 'cutoff'=500, 'thresh'=.04, 'Q'=100, 'iters'=4))
pars <- parlist[[5]]

f1 <- function(setNumber, cf.var, sc.var, nodes, mix, real.dat, run, pars, X){
  dpath <- paste('parOpt/dat/data-', setNumber, sep='')
  rpath <- paste('parOpt/res/res-', setNumber, sep='')
  cl <- makeForkCluster(nodes)
  registerDoParallel(cl)
  foreach(x=unlist(X),.export=c('mix','real.dat', 'kmax', 'psisets','etasets','cnmodels', 'runAlgs', 'convert.exp','setNumber','parlist',
                                'psiOptim2','psiFilter','filter','dpath','rpath','filter.mut','runAlg','ddirichlet',
                                'rdirichlet','dbeta2','estBetaParams'), .packages=c('gtools')) %dopar% {
    #mu <- metadata$data.params$mu
    runAlgs(x, set=setNumber, cn=TRUE, minim.vaf=.04, datapath=dpath, respath=rpath, mu=70, cf=cf.var, 
            sc=sc.var, exp=FALSE, mix=mix, real.dat=real.dat, run=run, pars)                      
  }
  parLapply(cl,6:33,function(x){
    runAlgs(x, set=setNumber, cn=TRUE, minim.vaf=.04, datapath=dpath, respath=rpath, mu=70, cf=cf.var, 
            sc=sc.var, exp=FALSE, mix=mix, real.dat=real.dat, run=run, pars)
  })
}

f1('CLL', cf.var=TRUE, sc.var=FALSE, nodes=28, X=list(1:262), mix=FALSE, real.dat=TRUE, run=5, pars=parlist[[5]])
#Best pars so far:
#list('sigma0'=.4, 'theta'=.3, 'ktheta'=.5, 'ntheta'=.5, 'stheta'=.5, 'mtheta'=.5, 'alpha'=.5, 'thresh'=0)
#list('sigma0'=.4, 'theta'=.5, 'ktheta'=.5, 'ntheta'=.5, 'stheta'=.5, 'mtheta'=.5, 'alpha'=.5, 'thresh'=0)
#Maybe: list('sigma0'=.4, 'theta'=.7, 'ktheta'=.9, 'ntheta'=.5, 'stheta'=.5, 'mtheta'=.5, 'alpha'=.5, 'thresh'=0

#Mixtures:
f1('mix', cf.var=TRUE, sc.var=FALSE, nodes=28, X=list(1:300), mix=TRUE, real.dat=FALSE, run=NULL, pars=parlist[[5]])


###For parameter optimization of algorithm.
source('parOpt/clonefinder/algorithm.R')
sampleIndices <- t(sapply(1:10,function(i){
  seq(from=i,to=300,by=12)
}))
sigma.lrr.range <- c(.01,.2)
sigma.baf.range <- c(.01,.1)
theta.range <- c(.01,.99)
ktheta.range <- c(.01,.99)
ranges <- c(sigma.lrr.range[2]-sigma.lrr.range[1],sigma.baf.range[2]-sigma.baf.range[1],
            theta.range[2]-theta.range[1],ktheta.range[2]-ktheta.range[1])
pars1 <- list('sigma.lrr'=sigma.lrr.range[1],'sigma.baf'=sigma.baf.range[1],'theta'=theta.range[1],
           'ktheta'=ktheta.range[1],'mtheta'=.9, 'alpha'=.5, 'thresh'=.04, 'cutoff'=100, 'Q'=100, 'iters'=4)
iterations.1 <- 4
iterations.2 <- 3
dpath <- 'parOpt/dat/testdata/'

binOpt <- function(i){
  cl <- makeForkCluster(length(sampleIndices[[i]]))
  registerDoParallel(cl)
  output <- foreach(x=unlist(sampleIndices[i,]),.export=c('kmax', 'psisets','etasets','cnmodels', 'runAlgs',
      'psiOptim2','psiFilter','filter','dpath','filter.mut','runAlg','ddirichlet',
      'rdirichlet','dbeta2','estBetaParams','pars1','ranges','iterations.1','iterations.2'), .packages=c('gtools')) %dopar% {
        dat <- get(load(paste('parOpt/dat/testdata/mixdat-',x,'.rda', sep='')))
        sim <- get(load(paste('parOpt/sim/testmixtures/mixsim-',x,'.rda', sep='')))
        truePsi <- sim$psi
        trueK <- length(which(truePsi>0))
        parset1 <- pars1
        cndata <- dat$cn.data
        res.cf.base <- runAlg(cndata, vardata=NULL, cnmodels, psisets[[2]], pars=parset1, priorbank=NULL)
        kdist.base <- abs(ks-trueK)
        psidist.base <- sum((res.cf$psi-truePsi)^2)^.5
        for(z in 1:iterations.2){
          for(n in 1:4){
            interval <- ranges[n]
            for(j in 1:iterations.1){
              newpars <- parset1
              newpars[n] <- parset1[n] + interval
              res.cf <- runAlg(cndata, vardata=NULL, cnmodels, psisets[[2]], pars=newpars, priorbank=NULL)
              ks <- length(which(res.cf$psi>0))
              kdist <- abs(ks-trueK)
              psidist <- sum((res.cf$psi-truePsi)^2)^.5
              if(psidist<psidist.base){
                parset1 <- newpars
                sign <- -1*sign
              }
              interval <- sign*ranges[n]/j
            }
            parset1 <- newpars
          }
        }
        list('pars'=pars,'kdist'=kdist,'psidist'=psidist)
      }
  output
}

temp <- binOpt(1)
save(temp,file='parOpt/binoptim.rda')

###Uniform sampling:
constants <- c('mtheta'=.9, 'alpha'=.5, 'thresh'=.04, 'cutoff'=100, 'Q'=100, 'iters'=4)
g <- function(index.list,j){
  sigmas.lrr <- runif(j,sigma.lrr.range[1],sigma.lrr.range[2])
  sigmas.baf <- runif(j,sigma.baf.range[1],sigma.baf.range[2])
  thetas <- runif(j,theta.range[1],theta.range[2])
  kthetas <- runif(j,ktheta.range[1],ktheta.range[2])
  pardf <- data.frame('sigma.lrr'=sigmas.lrr,'sigma.baf'=sigmas.baf,'theta'=thetas,'ktheta'=kthetas)
  par.index <- rep(1:j,length(index.list))
  cl <- makeForkCluster(length(index.list))
  output <- parLapply(cl,1:(j*length(index.list)),function(k){
    x <- ceiling(k/length(i))
    y <- par.index[k]
    pars <- pardf[y,]
    pars <- as.list(c(pars,constants))
    dat <- get(load(paste('parOpt/dat/testdata/mixdat-',x,'.rda', sep='')))
    cndata <- dat$cn.data
    sim <- get(load(paste('parOpt/sim/testmixtures/mixsim-',x,'.rda', sep='')))
    truePsi <- sim$psi
    trueK <- length(which(truePsi>0))
    pars <- c(constants,pardf[y,])
    res.cf <- runAlg(cndata, vardata=NULL, cnmodels, psisets[[2]], pars=pars, priorbank=NULL)
    kdiff <- length(which(res.cf$psi>0))-trueK
    psi <- res.cf$psi
    if(length(psi)>length(truePsi)){
      truePsi <- c(truePsi,rep(0,length(psi)-length(truePsi)))
    }
    psidist <- sum((psi-truePsi)^2)^.5
    list('pars'=pars,'k'=length(which(res.cf$psi>0)),'kdiff'=kdiff,'psidist'=psidist,'x'=x,'y'=y)
  })
  output
}

t1 <- Sys.time()
run1 <- g(index.list=sampleIndices[1,],20)
t2 <- Sys.time()
rt <- difftime(t2,t1)
save(run1,file='parOpt/run1.rda')

