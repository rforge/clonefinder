#source("C:\\Users\\Mark\\Documents\\Lab\\Code\\Code-10-30\\objs4.R")
source("E:\\Current Files 5\\objs4.R")
nSeg <- 1000
wts <- rev(3^(1:5))
wts <- wts/sum(wts)
xy <- data.frame(x = c(.2, .7, .8, .1, .4),
                 y = c(.2, .3, .5, .9, .7))
library(mc2d)

# generate the markers explicitly
markers <- round(runif(nSeg, 25, 1000))
fracs<-c(4,3,1)
fracs<- sort(fracs, decreasing=TRUE)
TrueNclones<- length(fracs)
#nclones is the number of clones we beleive there to be
nclones<- 3
s0<- .1
truepsis<- fracs/sum(fracs)

# now simulate a tumor; length of 'fracs' in first argument is number of clones
abstractTumor <- AbstractTumor(fracs, markers, wts)
# and get the concrete representation
tumor <- Tumor(abstractTumor, xy)

logLikely <- function(dataset, phi, object, sigma0=.25) {
  # dataset = matrix with 'x' and 'y' columns produced by "generateData"
  # phi = vector of probabilities for each compartment
  # object = TumorByCompartment
  xy <- object@pureCenters
  markers <- object@markers
  if (length(phi) < nrow(xy)-1) stop("You did not supply enough entries in 'phi'")
  if (length(phi) > nrow(xy)) stop("You supplied too many entries in 'phi'")
  if (length(phi) < nrow(xy)) {
    lastphi <- 1 - sum(phi)
    phi <- c(phi, lastphi)
  }
  if (any(phi < 0)) stop("Negative probabailities are not allowed")
  phi <- matrix(phi/sum(phi), nrow=1) # make sure they add up to 1
  center <- as.data.frame(phi %*% as.matrix(xy))
  secondMomentX <- sum(phi*(xy[,1]^2+sigma0^2))
  secondMomentY <- sum(phi*(xy[,2]^2+sigma0^2))
  sigmaX <- (secondMomentX - (sum(phi*xy[,1]))^2)/sqrt(markers)
  sigmaY <- (secondMomentY - (sum(phi*xy[,2]))^2)/sqrt(markers)
  px <- dnorm(dataset$x, center$x, sigmaX, log=TRUE)
  py <- dnorm(dataset$y, center$y, sigmaY, log=TRUE)
  px+py
}

nPhi <- 10000

foo <- function(seed){
  set.seed(seed)
  simdata <- generateData(tumor, s0)
  
  phiset <- sampleSimplex(nPhi, 5)
  temp <- phiset
  
  temp[,1] <- temp[,1]-1
  #Not sure I understand the above line.
  euclid <- apply(temp^2, 1, sum)
  phiset <- phiset[order(euclid),]
  
  #phiset2<- read.table(file="F:\\Current Files 5\\phiset2.R")
  #write.table(phiset, file="F:\\Current Files 5\\phiset.R")
  #phiset<- read.table(file="F:\\Current Files 5\\phiset.R")
  
  #Compute all likelihoods.
  likelihoods <- apply(phiset, 1, function(phi) {
    logLikely(simdata, phi, tumor, sigma0=s0)
  })
  maxLikeIndex <- apply(likelihoods, 1, which.max)
  phipick <- phiset[maxLikeIndex,]
  
  # pass2
  
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
  rm(i, index, iset, temp, euclid)
  ls()
  
  # get the likelihoods for the new phis
  likelihoods <- apply(newphiset, 1, function(phi) {
    logLikely(simdata, phi, tumor, sigma0=s0)
  })
  maxLikeIndex <- apply(likelihoods, 1, which.max)
  
  phipick <- newphiset[maxLikeIndex,]
  phiv <- as.vector(phipick)
  l<- 2^nclones
 
  denx<- density(phiv, bw="SJ")$x
  boole1 <- denx>=0
  boole2 <- denx<=1
  cutoff1<- length(which(boole1 %in% FALSE))
  cutoff2<- length(which(boole1 %in% FALSE))
  #ind<- which(boole %in% TRUE)
  den<- density(phiv, bw="SJ")$y
  #den<- density(phiv)$y[(cutoff1+1):(length(denx)-cutoff2-1)]
  plot(den)
  plot(density(phiv, bw="SJ"))
  
  smoothfun<- function(vec){
    innerfun<- function(i){
      mean(vec[(i-5):(i+5)])
    }
    sapply(6:(length(vec)-5), innerfun)
  }
  smoothed<- smoothfun(den)
  plot(smoothed)
  den <- smoothed
  
  densfun<-function(vec){
    ifun<-function(i){
      vec[i]>vec[i-1] & vec[i]>vec[i+1]
    }
    unlist(sapply(1:length(den), ifun))
  }
  indices<- unlist(densfun(den))
  maxima<- (which(indices %in% TRUE)+1)
  maxima<- (maxima[2:(length(maxima)-1)]-cutoff1)/(length(denx)-cutoff1-cutoff2)
  densities<- den[(which(indices %in% TRUE)+1)][2:(length(maxima)+1)]
  
  list(data.frame(maxima, densities), phipick)
}
dd <- foo(1)
df<- dd[[1]]
phipick <- dd[[2]]
j <- 2^nclones
peaks<-df[with(df, order(densities)), ][(nrow(df)-j):nrow(df),]
#note: '-7' means I took the maxima with the top 8 densities.
#There should only be 6, but I added 2 just in case there 
#are duplicates of one or two maxima. 2 is more or less an
#arbitrary choice.This was because it is better to retain 
#arbitrary data than lose necessary data.
#peaks<- peaks[peaks$densities>1,]
peaks<-peaks[with(peaks, order(maxima)), ]

###clustering and stuff
km <- kmeans(df$maxima,2^nclones+1)
clusters <- km$cluster
clustdf <- cbind(df, clusters)
wtfun <- function(c){
  d<- clustdf[clustdf$clusters==c,]
  weights<- d$densities/sum(d$densities)
  sum(weights*d$maxima)
}
weightedmeans <- sort(sapply(1:max(clusters), wtfun))
phis=sapply(1:length(truepsis), function(n){combn(truepsis, n, sum)})
phis=sort(unlist(phis))

#vec <- peaks$maxima
vec <- weightedmeans
indices <- combn(1:(length(vec)), nclones)
#indices <- expand.grid(rep(list(1:length(vec)), nclones)) 
#indices <- t(as.matrix(indices))
innerfun <- function(n){
    psis <- vec[unlist(indices[,n])]
    combos <- lapply(2:(nclones-1), function(x){combn(psis, x)})
    set <- vec[-indices[,n]]
    indices2 <- combn(1:length(set), ncol(combos[[1]]))
    #indices2 <- expand.grid(rep(list(1:length(set)), ncol(combos[[1]]))) 
    #indices2 <- t(as.matrix(indices2))
    matchfun <- function(i){
    sum((set[indices2[,i]]-colSums(combos[[1]]))^2)
    }
    val<- min(sapply(1:(ncol(indices2)), matchfun))
    if(nclones==2){
      val<- 0
    }
    val + (1-sum(psis))^2
}
q <- which.min(sapply(1:ncol(indices), innerfun))
ivec <- c(indices[,q])
psis <- sort(vec[ivec])

data.frame(psis=rev(psis), truepsis)
data.frame(psis=sort(rev(psis)/sum(rev(psis)), decreasing=TRUE), truepsis)

#data.frame(means=sort(weightedmeans), phis[-(2^nclones-1)])
#phis
#weightedmeans


#Started to 
phiv <- as.vector(phipick)
vec1 <- sapply(1:5000, function(i){which.min((phiv[i]-combos)^2)})

foo2 <- function(clone){
matfun<- function(n){
  k <- vec1[n]
  v <- binary[k]
  str<- toString(v)
  cvec <- strsplit(str, "")
  cvec<- rev(as.numeric(cvec[[1]]))
  w <- which(cvec %in% 1)
  cvec[clone]==1
}
t(sapply(1:5000, matfun))
}
boole <- c(unlist(foo2(1)))
boole[boole==TRUE]<- 1
boole[boole==FALSE]<- 0
zed1 <- matrix(boole, nrow=1000, ncol=5)
