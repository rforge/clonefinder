#################################################
source("objs4.R")
set.seed(363453)

# Parameters that we need to define the structure
nSeg <- 1000
wts <- rev(5^(1:5))
wts <- wts/sum(wts)
xy <- data.frame(x = c(.2, .7, .8, .1, .4),
                 y = c(.2, .3, .5, .9, .7))

# generate the markers explicitly
markers <- round(runif(nSeg, 25, 1000))
# now simulate a tumor
abstractTumor <- AbstractTumor(c(3,1), markers, wts)
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
# These aare not used by the model-fitting code.
Z1 <- makeLatent(tumor@data[,1])
Z2 <- makeLatent(tumor@data[,2])
PhiMatrix <- tumor@fraction[1]*Z1 + tumor@fraction[2]*Z2

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
phiv <- as.vector(phipick)
km <- kmeans(phiv, centers=c(0, 1/3, 2/3, 1))
plot(phiv, col=km$cluster, pch=16)
table(km$cluster)

# assuming each cluster is a separate beta distribution, estimate parametera
# by the method of moments
params <- as.matrix(sapply(1:4, function(i) estBetaParams(phiv[km$cluster == i])))
params
# the middle two should be symmetric, so let's force it
# (could proabbaly have been msrter about this earlier)
psis <- (params[1,] / apply(params, 2, sum))[2:3]
psi1 <- mean(c(psis[2], 1-psis[1]))
psi2 <- 1 - psi1
c(psi1, psi2)

# Now we can get the first estimate of the latent-variable z-matrices
z1map <- c(0, 0, 1, 1) # first clone is the 75% one
z2map <- c(0, 1, 0, 1) # second clone is the 25% one
mkm <- matrix(km$cluster, ncol=5)
zed1 <- matrix(z1map[km$cluster], ncol=5)
zed2 <- matrix(z2map[km$cluster], ncol=5)
# And we have a re-estimated phi matrix
phiest <- psi1*zed1 + psi2*zed2

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

mean(Z1 != zed1)
mean(Z2 != zed2)

# Get set up to perform EM algorithm to improve the estimates.
# Since psi ranges between 0 and 1, we are going to use a
# logistic transformation to work on the real line and avoid
# problems that might crop up near te boundary.
lp <- function(p) log(p/(1-p)) 
ea <- function(a) {
  temp <- exp(a)
  temp/(1+temp)
}

# Computes the NEGATIVE log-likelihood for some value of psi,
# assuming that we already have the z1 and z2 matrices.
# Note that the input is assumed to be given by
#      x = lp(psi)
myTarget <- function(x, z1, z2, data, tumor) {
  psi <- ea(x)
  phinew <- psi*z1 + (1-psi)*z2
  loglikes <- sum(tock <- sapply(1:nrow(phinew), function(i, phi) {
    sum(log(likely(data[i,], phi[i,], tumor, sigma0=0.25)))
  }, phi=phinew))
  - loglikes
}


# parameters that control the EM loop
currlike <- 0
lastlike <- -10^5
epsilon <- 100 # only small compared to thesize ofthe likelihood

while(abs(lastlike - currlike) > epsilon) {
# M-step: Given Z1 and Z2, use MLE to find optimal psi
  runner <- optimize(myTarget, c(-5,5), z1=zed1, z2=zed2, data=simdata, tumor=tumor)
  psi <- ea(runner$minimum)
  lastlike <- currlike
  currlike <- -runner$objective
  cat("Log likelihood: ", currlike, "\n", file=stderr())
# E-step: Given psi, get values for Z1 and Z2
  sigma0 <- 0.25
  holdme <- matrix(NA, ncol=25, nrow <- 1000) 
  # rows = segments. For columns,  need to consider 5*5 = 25 possible
  # assignments of pure compartments to segments.
  for (i in 1:5) {
    for (j in 1:5) {
      # given compartment assignments i and j, compute the weighted center
      xy <- psi*tumor@pureCenters[i,] + (1-psi)*tumor@pureCenters[j,]
      # then get the independent likelihoods
      dx <- dnorm(simdata$x, xy$x, sigma0/sqrt(tumor@markers))
      dy <- dnorm(simdata$y, xy$y, sigma0/sqrt(tumor@markers))
      # take the product, but use a log transform
      pp <- log(dx)+log(dy)
      holdme[, 5*(i-1) + j] <- pp
    }
  }
  # find the best assignments
  picker <- apply(holdme, 1, which.max)
  # next few lines have to back compute i and j from the column
  # number in 1..25.
  J <- picker %% 5
  J[J==0] <- 5
  I <- 1 + (picker - J)/5
  zed1 <- makeLatent(I)
  zed2 <- makeLatent(J)
  if (any(apply(zed1, 1, sum) != 1)) stop("bad zed1")
  if (any(apply(zed2, 1, sum) != 1)) stop("bad zed1")
}

# examine the results. We get almost everything correct
mean(Z1 == zed1)
mean(Z2 == zed2)
sum(Z1 != zed1)
sum(Z2 != zed2)
psi

# trying to find some way to visualize the results
foo <- (Z1 != zed1) | (Z2 != zed2)
cfoo <- apply(foo, 1, sum)

chunk <- 2*(zed2 > 0) + 1*(zed1>0)
oc <- order(chunk[,1], chunk[,2], chunk[,3], chunk[,4])

#oc <- order(tumor@compartments[,1], tumor@compartments[,2], tumor@compartments[,3], tumor@compartments[,4])
image(1:1000, 1:5, chunk[oc,], col=c("gray", "blue", "red", "magenta"),
      xlab="Segments", ylab="Compartments")

# add an extra row that flagsd deviations from the truth.
# can make N relatively small to focus on the errors.
N <- 1000
image(1:N, 1:6, cbind(tumor@compartments, 0.5*(cfoo>0))[oc,][1:N,], col=c("gray", "blue", "yellow", "red", "magenta"),
      xlab="Segments", ylab="Compartments")
abline(h=5.5)
