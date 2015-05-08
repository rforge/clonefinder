#################################################
source("C:\\Users\\Mark\\Documents\\Lab\\Code\\Code-8-15\\objs4.R")
#source("F:\\Current Files 3\\objs4.R")
#source("objs4.R")
set.seed(363453)

# Parameters that we need to define the structure
TrueNclones<- 3
nSeg <- 1000
wts <- rev(5^(1:5))
wts <- wts/sum(wts)
xy <- data.frame(x = c(.2, .7, .8, .1, .4),
                 y = c(.2, .3, .5, .9, .7))

# generate the markers explicitly
markers <- round(runif(nSeg, 25, 1000))
fracs<-c(3,2,1)
TrueNclones<- length(fracs)
#nclones is the number of clones we beleive there to be
nclones<- 3

# now simulate a tumor; length of 'fracs' in first argument is number of clones
abstractTumor <- AbstractTumor(fracs, markers, wts)
# and get the concrete representation
tumor <- Tumor(abstractTumor, xy)
# clean up by removing stuff we don't need
rm(markers, nSeg, wts, xy, abstractTumor)
ls()

# simulate data by selecting the weighted means with appropriate
# standard error of the mean
simdata <- generateData(tumor)

# example: compute likelihood of equal weight of compartments
phi <- rep(0.2, 4)
temp <- likely(simdata, phi, tumor)
summary(temp)
rm(phi, temp)

# simulate a uniform set of phi-vectors
nPhi <- 10000
phiset <- sampleSimplex(nPhi, 5)
# check that it makes sense
summary(apply(phiset, 1, sum))
opar <- par(mfrow=c(2,1))
plot(density(phiset[, 5]))
hist(phiset[,5], breaks=123)
par(opar)
rm(opar)
# Reorder the rows by distance from (1,0,...,0).
# This is only wiorth doing since we know that the most
# common correct answer should be closest to here.
temp <- phiset
temp[,1] <- temp[,1]-1
euclid <- apply(temp^2, 1, sum)
phiset <- phiset[order(euclid),]
rm(temp, euclid)

# compute all the likelihoods
likelihoods <- apply(phiset, 1, function(phi) {
  likely(simdata, phi, tumor, sigma0=0.25)
})
class(likelihoods)
dim(likelihoods)

maxLikeIndex <- apply(likelihoods, 1, which.max)
plot(maxLikeIndex)
# Note: as expected, most segments are of the form [1,1], so
# both clones use the same first compartment. So, most times
# we should prefer the phi-set closest to (1, 0, ..., 0)

phipick <- phiset[maxLikeIndex,]
dim(phipick)
hist(phipick, breaks=123)

#########################
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
rm(i, index, iset, temp, euclid)
ls()

# get the likelihoods for the new phis
likelihoods <- apply(newphiset, 1, function(phi) {
  likely(simdata, phi, tumor, sigma0=0.25)
})
maxLikeIndex <- apply(likelihoods, 1, which.max)
plot(maxLikeIndex)

phipick <- newphiset[maxLikeIndex,]
hist(phipick, breaks=123)
# Now we see bumps near 0.25 and 0.75, as we should

####################################################
# next step is to find a formal way to "discover" the
# bumps that we can see visually. 

# First, we record the truth about the segments.
# These are not used by the model-fitting code.

zfun<-function(n){
  makeLatent(tumor@data[,n])
}
ZMatrices<-lapply(1:TrueNclones, zfun)
ZMatrices<-rev(ZMatrices)

PhiMatFun<-function(n){
  tumor@fraction[n]*ZMatrices[[n]]
}
PhiMatrix <- Reduce('+', lapply(1:TrueNclones, PhiMatFun))

# model fitting
# idea is to fit the set of (esimtated) phi-values as a
# mixture of four beta distributions
#     Beta(a, 1), a > 1. Peaked near 1
#     Beta(1, b), b > 1. Peaked near 0
#     Beta(c, d), c > 1, d > !.  Peaked in between
#     Beta(d, c), Symmetric peak
if(FALSE) { # plots to show that we have the right kind of betas
  xx <- seq(0, 1, length=501)
  yy <- dbeta(xx, 10, 1)
  plot(xx, yy,type='l')

  yy <- dbeta(xx, 1, 40)
  plot(xx, yy,type='l')
  rm(xx, yy)
}

# start by using k-means to approximate the four groups
#Discerning number of clones: within cluster sum of squares threshold?
centervec<-c(seq(0, 1, length=nclones+2))
phiv <- as.vector(phipick)
km <- kmeans(phiv, centers=centervec)
plot(phiv, col=km$cluster, pch=16)
table(km$cluster)
km[2]

# assuming each cluster is a separate beta distribution, estimate parameters
# by the method of moments
params <- as.matrix(sapply(1:(nclones+2), function(i) estBetaParams(phiv[km$cluster == i])))
params

#Replaced lines  that are peculiar to 2 clones situation; psi1 is psis[1] for example
psis <- (params[1,] / apply(params, 2, sum))[2:(nclones+1)]
psis<-rev(psis)
symmfun<-function(n){
  mean(c(psis[n], 1-sum(psis[-n])))
}
psis<-sapply(1:nclones, symmfun)
# Now we can get the first estimate of the latent-variable z-matrices

zmapfun<-function(clone){
zmap<-c(rep(0, nclones+1), 1)
zmap[nclones+1-(clone-1)]<-1
zmap
}
#arithmetic explanation: there are nclones+2 clusters; if clones are in increasing order, then the '1' for 
#clone 1 in zmap should be at nclones+2-1-(clone number-1)=nclones+1-0; or nclones+1-(clone-1)
zmapfun(2)
zmapList<-sapply(1:nclones, zmapfun)
#so now, for example, z2map is zmapList[,2]

mkm <- matrix(km$cluster, ncol=5)

zedfun<-function(clone){
  matrix(zmapList[,clone][km$cluster], ncol=5)
}
zedlist<-lapply(1:nclones, zedfun)
#now zed1 is zedlist[[1]]

# And we have a re-estimated phi matrix
#phiest <- psi1*zed1 + psi2*zed2

phiestfun<-function(clone){
  psis[clone]*zedlist[[clone]]
}

phiestmatrix<-lapply(1:nclones, phiestfun)

#Get a 'summed list' by summing each element in the list of lists
phiest<-Reduce("+", phiestmatrix)

# how close was the old estimate to the truth?
summary(PhiMatrix - phipick)
hist(PhiMatrix - phipick, breaks=57)
mean(abs(PhiMatrix - phipick) > 0.05) # got almost 94% of them right

# How much did the estimate change?
summary(phipick - phiest)
hist(phipick - phiest, breaks=57)
mean(abs(phipick - phiest) > 0.05)

# how close is the new estimate to the truth?
summary(PhiMatrix - phiest)
hist(PhiMatrix - phiest, breaks=57)
image(PhiMatrix - phiest)
mean(abs(PhiMatrix - phiest) > 0.05) # got almost 99% of them right

# Get set up to perform EM algorithm to improve the estimates.
# Since psi ranges between 0 and 1, we are going to use a
# logistic transformation to work on the real line and avoid
# problems that might crop up near the boundary.
lp <- function(p) log(p/(1-p)) 
ea <- function(a) {
  temp <- exp(a)
  temp/(1+temp)
}

library(Brobdingnag)

# Computes the NEGATIVE log-likelihood for some value of psi,
# assuming that we already have the z1 and z2 matrices.
# Note that the input is assumed to be given by
#      x = lp(psi)
#Current issue: 'ea' has been removed, see what happens...
myTarget <- function(par, zlist, data, tumor) {
  par<-sort(par, decreasing=TRUE)
  psiset <- ea(par)
  productfun<-function(n){
    psiset[n]*zlist[[n]]
  }
  productlist<-lapply(1:(nclones-1), productfun)
  phinew<-Reduce("+", productlist)+(1-sum(psis[1:(nclones-1)]))*zlist[[nclones]]
  loglikes <- sum(tock <- sapply(1:nrow(phinew), function(i, phi) {
    sum(logLikely(data[i,], phi[i,], tumor, sigma0=0.25))
  }, phi=phinew))
  - loglikes
}
myTarget(par=c(.5, .333), zedlist, simdata, tumor)

#source("C:\\Users\\Mark\\Documents\\Lab\\Code\\Code-9-25\\Z.R")
#mean(zed2==zedlist[[2]])

# parameters that control the EM loop
currlike <- 0
lastlike <- -10^5
epsilon <- 100 # only small compared to the size ofthe likelihood

library(gtools)

while(abs(lastlike - currlike) > epsilon) {
# M-step: Given Z1 and Z2, use MLE to find optimal psi
  startvec<-c(rep(1/nclones, nclones-1))
  runner<-optim(par=startvec, myTarget, zlist=zedlist, data=simdata, tumor=tumor)
  psivec <- c(ea(c(runner$par[1:(nclones-1)])), 1-sum(ea(c(runner$par[1:(nclones-1)]))))
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
  holdme <-sapply(1:nrow(perms), rowfun, psivals=psis)
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
resultfun<-function(n){
  paste("Clone", n, ":", "mean of Z=zed:", mean(ZMatrices[[n]] == zedlist[[n]]), "; sum of false assignment:", 
        sum(ZMatrices[[n]] != zedlist[[n]]), sep=" ")
}
lapply(1:TrueNclones, resultfun)
psis/sum(psis)

##############Below: to be modified
# trying to find some way to visualize the results
foo <- (Z1 != zed1) | (Z2 != zed2)
cfoo <- apply(foo, 1, sum)

chunk <- 2*(zed2 > 0) + 1*(zed1>0)
oc <- order(chunk[,1], chunk[,2], chunk[,3], chunk[,4])

#oc <- order(tumor@compartments[,1], tumor@compartments[,2], tumor@compartments[,3], tumor@compartments[,4])
image(1:1000, 1:5, chunk[oc,], col=c("gray", "blue", "red", "magenta"),
      xlab="Segments", ylab="Compartments")

# add an extra row that flags deviations from the truth.
# can make N relatively small to focus on the errors.
N <- 1000
image(1:N, 1:6, cbind(tumor@compartments, 0.5*(cfoo>0))[oc,][1:N,], col=c("gray", "blue", "yellow", "red", "magenta"),
      xlab="Segments", ylab="Compartments")
abline(h=5.5)

