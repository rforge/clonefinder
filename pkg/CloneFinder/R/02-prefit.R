# first pass at fitting the model

setClass("PrefitCloneModel",
         representation=list(
           data="data.frame",
           phiset = "matrix",
           likelihoods = "matrix",
           phipick="matrix",
           maxLikeIndex = "vector"))

.reorderVectors <- function(phiset) {
# Reorder the rows by distance from (1,0,...,0).
# This is only worth doing since we know that the most
# common correct answer should be closest to here.
  temp <- phiset
  temp[,1] <- temp[,1] - 1
  euclid <- apply(temp^2, 1, sum)
  phiset[order(euclid),]
}

.computeLikelihoods <- function(phiset, segmentdata, tumor, log=TRUE) {
  apply(phiset, 1, function(phi) {
    likely(segmentdata, phi, tumor, log)
  })
}

PrefitCloneModel <- function(segmentdata, nPhi = 10000) {
# simulate a uniform set of phi-vectors
  phiset <- .reorderVectors(sampleSimplex(nPhi, 5))
  likelihoods <- .computeLikelihoods(phiset, segmentdata, tumor, TRUE)
# locate the maximum likelihood for each phi-vector
  maxLikeIndex <- apply(likelihoods, 1, which.max)
  phipick <- phiset[maxLikeIndex,]
# return the results, bundled in an S4 class
  new("PrefitCloneModel",
      data=segmentdata,
      phiset=phiset,
      likelihoods=likelihoods,
      maxLikeIndex = maxLikeIndex,
      phipick = phipick)
}

updatePhiVectors <- function(object, nPhi=10000) {
  if (!inherits(object, "PrefitCloneModel")) {
    object <- PrefitCloneModel(object, nPhi)
  }
  multiplier <- round(nPhi/1000)
# resample the phi vectors to be near the ones selected as
# optimal in the first pass
  newphiset <- matrix(NA, ncol=5, nrow=nPhi)
  for (i in 1:1000) {
    index <- 1 + multiplier*(i-1)
    iset <- index:(index+multiplier-1)
    newphiset[iset,] <- rdirichlet(multiplier, 2*object@phiset[object@maxLikeIndex[i],])
  }
  newphiset <- .reorderVectors(newphiset)
# get the likelihoods for the new phis
  likelihoods <- .computeLikelihoods(newphiset, object@data, tumor, TRUE)
  maxLikeIndex <- apply(likelihoods, 1, which.max)
  phipick <- newphiset[maxLikeIndex,]
  new("PrefitCloneModel",
      data=object@data,
      phiset=newphiset,
      likelihoods=likelihoods,
      maxLikeIndex = maxLikeIndex,
      phipick = phipick)
}

setMethod('plot', signature(x='PrefitCloneModel', y='missing'),
function(x, xlab="Set of Phi Vectors", ylab="Index to Maximum Likelihood", ...) {
  plot(x@maxLikeIndex, xlab=xlab, ylab=ylab, ...)
})

setMethod('hist', signature(x='PrefitCloneModel'),
function(x, xlab="Cell Fraction (Phi)", main="", ...) {
  if(is.null(main)) main <- deparse(substitute(x))
  hist(x@phipick, xlab=xlab, main=main, ...)
})

setMethod('summary', signature('PrefitCloneModel'), function(object, ...) {
  summary(apply(object@phiset, 1, sum))
})
