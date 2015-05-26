# Want to avoid a constrained optimization.  But the parameters of he model
# are given by a vector psi=(psi_1, .., psi_N) that satisifies \sum psi_i = 1
#
# so, we want to transform an unconstrained set of either N or preferably N-1
# independent variables.  Here the idea is to give preference to the first entries
# so that: 0 < x1 < 1,   0 < x2 < 1-x1,    0 < x3 < 1-x1-x2, etc.
# and rewrite  Vi = xi/(1 - (x1+x2+x{i-1}))
# to reverse, xi = v1(1-v{i-1})(...(1-v2)(1-v1)
#
# These are computed by the "forward" and "backward" transformations.

ea <- function(a) {z <- exp(a); z/(1+z)}
lp <- function(p) log(p/(1-p))

# given an N-tuple of nonnegative values that sum to 1, convert it to
# an (N-1)-tuple of arbitrary real numbers
forward <- function(F) {
  L <- ncol(F)-1
  bonk <- 1 - t(apply(F, 1, cumsum))
  b <- as.matrix(data.frame(1, bonk[, 1:L]))
  lp(F/b)[, 1:L]
}

# Given an L-tuple of arbitrary real numbers, convert it to an
# (L+1)-tuple of non-negative reals that sum to 1
backward <- function(w) {
  w <- ea(w)
  weird <- as.matrix(data.frame(1, t(apply(1-w, 1, cumprod))))
  as.matrix(data.frame(w,1))*weird
}

# check forward and backward
if(FALSE) {
  x <- runif(10000, 0, 1)
  V <- matrix(lp(x), ncol=5)
  F <- t(apply(V, 1, forward)) # nonneg and sum to 1
  temp <- apply(F, 1, sum)

  w <- forward(F)
  flap <- backward(w)
  summary(as.vector(F - flap))
}



# Computes the NEGATIVE log-likelihood for some value of psi,
# assuming that we already have the Z matrices.
# Note that the input is assumed to be given by
#      x = forward(psi)   (so psi=backward(x))
# and Zs = 3D array of z-matrices, so should have dimension
#     nSegments x nCompartments x nClones
#     where nClones equals length(x)
myTarget <- function(x, Zs, data, tumor) {
  if (length(x) != dim(Zs)[3]) stop('mismatched argmuent sizes')
  psi <- backward(x)
  # next formula is vectorization of " phi = sum_i (psi_i * Z_i) "
  phinew <- apply(sweep(Zs, 3, psi, "*"), 1:2, sum)
  loglikes <- sum(tock <- sapply(1:nrow(phinew), function(i, phi) {
    sum(log(likely(data[i,], phi[i,], tumor, sigma0=0.25)))
  }, phi=phinew))
  - loglikes  # negate
}

# check algorithm
if (FALSE) {
  Zs <- array(1:60, dim=c(4, 5, 3))
  psi <- 1:3
  y <- sweep(Zs, 3, psi, "*")
  phinew <- apply(y, 1:2, sum)
}


# precompute all possible compartment-clone assignments
precomputeZed <- function(nComp, nClone) {
  base <- matrix(0, nrow=nComp, ncol=nComp)
  diag(base) <- 1
  ary <- array(NA, dim=c(nComp^nClone, nComp, nClone))
  dex <- 1:nComp
  for (i in 1:nClone) {
    edex <- rep(rep(dex, each=nComp^(i-1)), times=nComp^(nClone-i))
    ary[,,i] <- base[edex,]
  }
  ary
}

#############################
# simulate a sample data set
library(CloneFinder)
set.seed(539121)
# pure centers
xy <- data.frame(x = log10(c(2, 2, 1, 3, 4)/2),
                 y = c(1/2, 0, 0, 1/3, 1/4))
# number of segments
nSeg <- 1000
# number of SNP markers per segment
markers <- round(runif(nSeg, 25, 1000))
compModel <- CompartmentModel(markers, xy, 0.25)
# probability of a pure cloncal segment in each compartment
wts <- rev(5^(1:5))
wts <- wts/sum(wts)
# percentage of cells in each (of three) clone(s)
psis <- c(0.6, 0.3, 0.1)

tumor <- Tumor(compModel, psis, wts)
dataset <- generateData(tumor)

# prefit the model
pcm <- PrefitCloneModel(dataset)
# update it
upd <- updatePhiVectors(pcm)
slotNames(upd)

nClones <- 3 # we know this because we're omniscient
zedary <- precomputeZed(5, nClones)

#########################
# need to make a not-stupid initial guess at the psi-values.
# the following is adapted from Mark's "Algorithm_New-v3".

library(combinat)
# pad a vector with trailing zeros to the desired length
# while loop is probably not very efficient, but who cares
pad <- function(v, d) {
  while(length(v) < d) v <- c(v, 0)
  v
}

# can view this as likelihood of estimated phi-values given psi
computeSSE <- function(psi, sphi) {
  ssfun <- function(x) sum(x^2)
  # check the three-way model
  res <- sweep(sphi, 2, psi, '-')
  ss <- apply(res, 1, ssfun) # per-segment sum-of-square errors
  # find all sums (of combinations) of psi's
  combos <- sort(unique(unlist(sapply(1:length(psi), function(i){combn(psi, i, sum)}))))
  # the rows of m are all the two-way models
  m <- t(sapply(1:length(psi), function(i){sort(c(psi[i], 1-psi[i]), decreasing=TRUE)}))
  padm <- t(apply(m, 1, pad, d=5)) # pad with zeros to number of compartments
  ss2 <- data.frame(ss, apply(padm, 1, function(mod2) {
                                res2 <- sweep(sphi, 2, mod2, '-')
                                apply(res2, 1, function(x) sum(x^2))
                              })) # SSE for two-way models
  sstotal <- sum(apply(ss2, 1, min)) # use the best model for each segment
  sstotal
}

guessPsi <- function(upd, nClones) {
  phiset <- upd@phipick # previously estimated percent of cells in each compartment
  # sort the phi-vectors with largest first
  sortedphi <- t(apply(phiset, 1, sort, decreasing=TRUE)) 
  # sample te space of possible psi-vectors ...
  field <- sampleSimplex(200, nClones)
  # ... and put them in decreasing order
  sfield <- t(apply(field, 1, sort, decreasing=TRUE))
  # pad them each out to the corect lenggth for the number of compartments
  padfield <- t(apply(sfield, 1, pad, d=5))
  # for each "padded field vector" of psi values, find the total SS for each segment
  sse <- apply(padfield, 1, computeSSE, sphi=sortedphi)
  pick <- which(sse == min(sse))
  sfield[pick,]
}

estpsi <- guessPsi(upd, 3)

setZs <- function(psi, zedary, simdata, tumor) {
  sigma0 <- 0.25
  holdme <- matrix(NA, ncol=nrow(zedary), nrow=nrow(simdata))
  # rows = segments. For columns,  need to consider 5*nClone possible
  # assignments of pure compartments to segments.
  for (dex in 1:dim(zedary)[1]) { 
    # given compartment assignments, compute the weighted center
    wts <- zedary[dex,,] %*% matrix(psi, ncol=1)
    xy <- as.data.frame(t(wts) %*% as.matrix(tumor@pureCenters))
    # then get the independent likelihoods
    dx <- dnorm(simdata$x, xy$x, sigma0/sqrt(tumor@markers))
    dy <- dnorm(simdata$y, xy$y, sigma0/sqrt(tumor@markers))
    # take the product, but use a log transform
    pp <- log(dx)+log(dy)
    holdme[, dex] <- pp
  }
  # find the best assignments
  picker <- apply(holdme, 1, which.max)
  zedary[picker, ,]
}

Zmats <- setZs(estpsi, zedary, dataset, tumor)

########################
# parameters that control the EM loop
currlike <- 0
lastlike <- -10^5
epsilon <- 100 # only small compared to the size of the likelihood

while(abs(lastlike - currlike) > epsilon) {
# M-step: Given Z-matrices, use MLE to find optimal psi
  runner <- optim(rep(0, dim(Zmats)[3]), myTarget, Zs=Zmats, data=simdata, tumor=tumor)
  psi <- backward(runner$minimum)
  lastlike <- currlike
  currlike <- -runner$objective
  cat("Log likelihood: ", currlike, "\n", file=stderr())
# E-step: Given psi, get values for Z-matrices
  Zmats <- setZs(psi, zedary, dataset, tumor)
  if (any(apply(Zmats, c(1,3), sum) != 1)) stop("bad Z")
}



