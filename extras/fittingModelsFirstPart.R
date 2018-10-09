source("E:\\Current Files 5\\objs4.R")
set.seed(363453)

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

simdata <- generateData(tumor, .25)

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
sum(logLikely(simdata, phiset[1,], tumor))

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
  logLikely(simdata, phi, tumor, sigma0=0.25)
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
  logLikely(simdata, phi, tumor, sigma0=0.25)
})
maxLikeIndex <- apply(likelihoods, 1, which.max)

phipick <- newphiset[maxLikeIndex,]
hist(phipick, breaks=123)

phiv <- as.vector(phipick)
centers<-c(seq(0, 1, length=nclones+2))
km <- kmeans(phiv, centers)
params <- as.matrix(sapply(1:(nclones+2), function(i) estBetaParams(phiv[km$cluster == i])))

psis <- (params[1,] / apply(params, 2, sum))[2:(nclones+1)]
psis<- psis/sum(psis)
truepsis<- rev(fracs/sum(fracs))

psis
truepsis

testphiset<- newphiset[1,]
sum(logLikely(simdata, testphiset, tumor))
########################

zmapfun<-function(clone){
  zmap<-c(rep(0, nclones+1), 1)
  zmap[nclones+1-(clone-1)]<-1
  zmap
}
zmapList<-sapply(1:nclones, zmapfun)

zedfun<-function(clone){
  matrix(zmapList[,clone][km$cluster], ncol=5)
}
zedlist<-lapply(1:nclones, zedfun)

mean(zedlist[[1]]==Z[[1]])
mean(zedlist[[2]]==Z[[2]])
mean(zedlist[[3]]==Z[[3]])