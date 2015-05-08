source("E:\\Current Files 5\\objs4.R")
set.seed(363453)

nSeg <- 1000
wts <- rev(5^(1:5))
wts <- wts/sum(wts)
xy <- data.frame(x = c(.2, .7, .8, .1, .4),
                 y = c(.2, .3, .5, .9, .7))

# generate the markers explicitly
markers <- round(runif(nSeg, 25, 1000))
fracs<-c(4,3,2,1)
fracs<- sort(fracs, decreasing=TRUE)
TrueNclones<- length(fracs)
#nclones is the number of clones we beleive there to be
nclones<- 4
s0<- .1

# now simulate a tumor; length of 'fracs' in first argument is number of clones
abstractTumor <- AbstractTumor(fracs, markers, wts)
# and get the concrete representation
tumor <- Tumor(abstractTumor, xy)
# clean up by removing stuff we don't need
ls()

#Truth:
ZFun <- function(n) {
  makeLatent(tumor@data[,n])
}
Z<- lapply(1:nclones, ZFun)
PhiMatrixFun <- function(n){
  tumor@fraction[n]*Z[[n]]
}
PhiMatrices <- lapply(1:nclones, PhiMatrixFun)
PhiMatrix<- Reduce("+", PhiMatrices)

simdata <- generateData(tumor, s0)

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
dim(phipick)
hist(phipick, breaks=123)

# pass2
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
rm(i, index, iset, temp, euclid)
ls()

# get the likelihoods for the new phis
likelihoods <- apply(newphiset, 1, function(phi) {
  logLikely(simdata, phi, tumor, sigma0=s0)
})
maxLikeIndex <- apply(likelihoods, 1, which.max)

phipick <- newphiset[maxLikeIndex,]
hist(phipick, breaks=123)

phiv <- as.vector(phipick)
l<- 2^nclones
centers<-c(seq(0, 1, length=l))
km <- kmeans(phiv, centers)
params <- as.matrix(sapply(1:(nclones+2), function(i) estBetaParams(phiv[km$cluster == i])))

psis <- (params[1,] / apply(params, 2, sum))[2:(nclones+1)]
psis<- psis/sum(psis)
truepsis<- rev(fracs/sum(fracs))
plot(phiv, col=km$cluster, pch=16)

psis
truepsis
prepsis<-psis

vars<- params[1,]*params[2,]/((params[1,]+params[2,])^2*(params[1,]+params[2,]+1))
vars<- vars[2:(nclones+1)]
sds<- vars^.5

sum(likelihoods[maxLikeIndex])
#back-computing:
#(a/(a+b))/sum(a/(a+b)) <- truepsis
#a and b are vectors

#Try computing psis another way
denx<- density(phiv)$x
boole1 <- denx>=0
boole2 <- denx<=1
cutoff1<- length(which(boole1 %in% FALSE))
cutoff2<- length(which(boole1 %in% FALSE))
#ind<- which(boole %in% TRUE)
den<- density(phiv)$y
#den<- density(phiv)$y[(cutoff1+1):(length(denx)-cutoff2-1)]
plot(den)
plot(density(phiv))
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

truepsis
maxima
densities

#rules
npeaks<- length(maxima)
if(npeaks==nclones-1){
  psis<- rep(maxima[1], nclones)
}else{
if(nclones==3 & npeaks==4){
  k<- which(densities %in% sort(densities, decreasing=TRUE)[1:2])
  psis<- c(rep(min(maxima[k]), 2), 1-2*min(maxima[k]))
}else{
  if(nclones==3 & npeaks>=5){
    psis<- c(maxima[1], maxima[2], 1-maxima[1]-maxima[2])
  }
}
}

#kmeans trial 2:
l<- 2^nclones
centers<-c(seq(0, 1, length=l))
km <- kmeans(phiv, centers)
params <- as.matrix(sapply(1:(nclones+2), function(i) estBetaParams(phiv[km$cluster == i])))
kmpsis <- (params[1,] / apply(params, 2, sum))[2:(nclones+1)]
kmpsis <- kmpsis/sum(kmpsis)
kmpsis
plot(phiv, col=km$cluster, pch=16)

library(combinat)
combfun<- function(m){
combn(psis, m, sum)
}
combos <- unlist(c(0, sapply(1:length(psis), combfun)))
#0, 1, 2, 3, 1+2, 1+3, 2+3, 1+2+3...)

vfun <- function(m){
  combn(psis, m)
}
vlist<- sapply(1:length(psis), vfun)

clustfun<-function(n){
which.min((combos-mean(phiv[which(km$cluster %in% n)]))^2)
}
indexlist<- sapply(1:l, clustfun)
Mat<- matrix(km$cluster, nrow=1000, ncol=5)

ZMatFun<-function(){
  ones <- which(Mat %in% v)
  m<- matrix(rep(0, 5000), nrow=1000m ncol=5)
  m[ones] <- 1
  m
}

#How to decompose the maxima into psis, given the number of clones, using maxima, densities, and wts?
########################
vec<- 10^(0:(nclones-1))
combofun2<-function(m){
combn(vec, m, sum)
}
binary<- unlist(c(0, sapply(1:(nclones), combofun2)))
binary <- sort(binary, decreasing=FALSE)

outerfun <- function(clone){
innerfun<-function(n){
  v <- binary[n]
  str<- toString(v)
  cvec <- strsplit(str, "")
  cvec<- rev(as.numeric(cvec[[1]]))
  w <- which(cvec %in% 1)
  cvec[clone]==1
}
boolean <- sapply(km$cluster, innerfun)
boolean[boolean==FALSE]<- 0
boolean[is.na(boolean)]<- 0
boolean[boolean==TRUE]<- 1
matrix(boolean, nrow=1000, ncol=5)
}
ZMatList <- lapply(1:nclones, outerfun)
#ZMatList[1][[1]] is zed1

resultfun <- function(clone){
  e<- (sum(ZMatList[clone][[1]]!=Z[clone][[1]])/2)/1000
  results <- paste("Rate of false assignment for clone ", clone, ": ", e, sep="")
  cat(results, "\n")
}
sapply(1:nclones, resultfun)

#Source of inaccuracy?

##########################
#Tasks (3, 1, 1): fit with 4, 6, 8, and decide which is best. See 'Silhouette width'; write overview 
#of algorithm; 

#########################################################

lp <- function(p) log(p/(1-p)) 
ea <- function(a) {
  temp <- exp(a)
  temp/(1+temp)
}

myTarget <- function(par, zlist, data, tumor) {
  par<-(sort(par, decreasing=FALSE))
  if(sum(par)>=1){
    par<- par/sum(par)
  }
  psiset <- ea(par)
  productfun<-function(n){
    psiset[n]*zlist[[n]]
  }
  productlist<-lapply(1:(nclones-1), productfun)
  phinew<- Reduce("+", productlist)+(1-ea(sum(par[1:(nclones-1)])))*zlist[[nclones]]
  loglikes <- sum(tock <- sapply(1:nrow(phinew), function(i, phi) {
    sum(logLikely(data[i,], phi[i,], tumor, sigma0=0.25))
  }, phi=phinew))
  - loglikes
}
myTarget(par=c(lp(.6), lp(.27)), ZMatList, simdata, tumor)

currlike <- 0
lastlike <- -10^5
epsilon <- 100 # only small compared to the size ofthe likelihood

library(gtools)
psivec <- psis[1:nclones-1]

while(abs(lastlike - currlike) > epsilon) {
  # M-step: Given Z1 and Z2, use MLE to find optimal psi
  # startvec<-c(rep(1/nclones, nclones-1))
  zedlist <- ZMatList
  runner<-optim(par=psivec, myTarget, zlist=zedlist, data=simdata, tumor=tumor, method="L-BFGS-B")
  psivec <- c(ea(c(runner$par[1:(nclones-1)])), 1-ea(sum(c(runner$par[1:(nclones-1)]))))
  psivec<- sort(psivec, decreasing=TRUE)
  #Stipulated that only if currlike is positive can lastlike take its value; for some reason I'm getting negative likelihoods.
  #if(currlike>0){
  lastlike <- currlike
  #}
  currlike <- -runner$value
  cat("Log likelihood: ", currlike, "\n", file=stderr())
  # E-step: Given psi, get values for Z1 and Z2
  sigma0 <- 0.25
  perms<- permutations(5, nclones, c(1:5), repeats.allowed=TRUE)
  rowfun<-function(row, psivals){
    xy <- matrix(psivals, nrow=1) %*% as.matrix(tumor@pureCenters[perms[row,],])
    xy <- as.data.frame(xy)
    dx <- dnorm(simdata$x, xy$x, sigma0/sqrt(tumor@markers))
    dy <- dnorm(simdata$y, xy$y, sigma0/sqrt(tumor@markers))
    # take the product, but use a log transform
    pp <- log(dx)+log(dy)
    pp
  }
  holdme <-sapply(1:nrow(perms), rowfun, psivals=psivec)
  # find the best assignments
  picker <- apply(holdme, 1, which.max)
  # next few lines have to back compute i and j from the column
  # number in 1..5^nclones.
  modulusfun<-function(clone){
    I<- ceiling(picker/(5^(clone-1))) %% 5
    I[I==0] <- 5
    I
  }
  cloneset<- lapply(1:nclones, modulusfun)
  
  zedfun2<-function(clone){
    makeLatent(cloneset[[clone]])
  }
  zedlist<-lapply(1:nclones, zedfun2)
}

# examine the results. We get almost everything correct
psivec/sum(psivec)
resultfun2<-function(n){
  paste("Clone", n, ":", "mean of Z=zed:", mean(Z[[n]] == zedlist[[n]]), "; sum of false assignment:", 
        sum(Z[[n]] != zedlist[[n]]), sep=" ")
}
lapply(1:nclones, resultfun2)