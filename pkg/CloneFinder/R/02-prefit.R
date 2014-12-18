# first pass at fitting the model

setClass("PrefitCloneModel",
         representation=list(
           phiset = "matrix",
           likelihoods = "matrix",
           phipick="matrix",
           maxLikeIndex = "vector"))

PrefitCloneModel <- function(segmentdata, nPhi = 10000) {
# simulate a uniform set of phi-vectors
  phiset <- sampleSimplex(nPhi, 5)
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
# locate the maximum likelihood for each phi-vector
  maxLikeIndex <- apply(likelihoods, 1, which.max)
  phipick <- phiset[maxLikeIndex,]
# return the results, bundled in an S4 class
  new("PrefitCloneModel", phiset=phiset,
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
