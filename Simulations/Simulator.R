source("E:\\Current Files 5\\objs4.R")
nSeg <- 1000
wts <- rev(5^(1:5))
wts <- wts/sum(wts)
xy <- data.frame(x = c(.2, .7, .8, .1, .4), y = c(.2, .3, .5, .9, .7))
markers <- round(runif(nSeg, 25, 1000))
#N1=number of sets of psi values generated per number of clones
N1 <- 10
#N2=number of replicate tumors generated per set of psi values
N2 <- 10
#Here, we have ten tumors per set, 10 sets per number of clones, four different
#numbers of clones (2, 3, 4, 5), so a total of 400 simulations generated.

psigen <- function(n){
#pool determines what values are possible for psis
  pool <- seq(.1, .9, length=17)
  samplegen <- function(seed){
    set.seed(seed)
    k <- sample(pool, n, replace=FALSE)
    sort(k, decreasing=FALSE)
  }
  lapply(1:N1, samplegen)
}
psivecs <- unlist(lapply(2:5, psigen), recursive=FALSE)

#SimGen is the function that actually generates the lists of objects:
SimGen <- function(i, data){
  fracs <- data[[i]]
  L <- length(fracs)
  truepsis <- sort(fracs, decreasing=FALSE)/sum(fracs) 
  innerfun <- function(seed){
    set.seed(seed)
  abstractTumor <- AbstractTumor(fracs, markers, wts)
  tumor <- Tumor(abstractTumor, xy)
  simdata <- generateData(tumor)
  zfun<-function(n){
    makeLatent(tumor@data[,n])
  }
  ZMatrices<-lapply(1:length(fracs), zfun)
 
  #Now, generating the phi estimates:
  nPhi <- 10000
  phiset <- sampleSimplex(nPhi, 5)
  temp <- phiset
  temp[,1] <- temp[,1]-1
  euclid <- apply(temp^2, 1, sum)
  phiset <- phiset[order(euclid),]
  likelihoods <- apply(phiset, 1, function(phi) {
    likely(simdata, phi, tumor, sigma0=0.25)
  })
  maxLikeIndex <- apply(likelihoods, 1, which.max)
  phipick <- phiset[maxLikeIndex,]
  
  # pass2
  # want to refine the estimates of the phis.
  library(mc2d)
  multiplier <- round(nPhi/1000)
  # resample the phi vectors to be near the ones selected as
  # optimal inthe first pass
  newphiset <- matrix(NA, ncol=5, nrow=nPhi)
  for (i in 1:1000) {
    index <- 1 + multiplier*(i-1)
    iset <- index:(index+multiplier-1)
    newphiset[iset,] <- rdirichlet(multiplier, 2*phiset[maxLikeIndex[i],])
  }
  # reorder them, for same reason as before
  temp <- newphiset
  temp[,1] <- temp[,1]-1
  euclid <- apply(temp^2, 1, sum)
  newphiset <- newphiset[order(euclid),]
  likelihoods <- apply(newphiset, 1, function(phi) {
   likely(simdata, phi, tumor, sigma0=0.25)
  })
  maxLikeIndex <- apply(likelihoods, 1, which.max)
  phipick <- newphiset[maxLikeIndex,]
  
  list("True psis"=truepsis,"Simulated data"=simdata,"Estimated phis"= phipick,
       "Z Matrices"= ZMatrices)
  }
  tumorset <- lapply(((N1+1):(N1+N2)), innerfun)
  filegen <- function(j){
    assign(paste("sim.", L, ".", i, ".", j, sep=""),tumorset[[j]])
    save(list=paste("sim.", L, ".", i, ".", j, sep=""), file=paste("F:\\Mark\\Simulations\\sim.",L, 
                                                                   ".", i,".",j,".rda", sep=""))
  }
  sapply(1:N2, filegen)
}
sapply(1:length(psivecs), SimGen, data=psivecs)
#Each .rda file contains a list with 4 objects: [1] true psi values, [2] simulated 'LRR' and 'BAF',
#[3] estimated phi values, and [4] true Z matrices, in order from smallest clone to largest.