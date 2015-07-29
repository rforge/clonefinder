# first pass at fitting the model

# KRC: Do we actually need all of the pieces we carry along here?
# KRC: Should we carry along the compartment model to reuse when we update phi vectors?
setClass("PrefitCloneModel",
         representation=list(
           data="data.frame",
           phiset = "matrix",
           likelihoods = "matrix",
           phipick="matrix",
           phiv="numeric",
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

.computeLikelihoods <- function(phiset, segmentdata, compartments, log=TRUE) {
  apply(phiset, 1, function(phi) {
    likely(segmentdata, phi, compartments, log)
  })
}

PrefitCloneModel <- function(segmentdata, compartments, nPhi = 10000) {
# recalibrate nPhi based on actual number of segments
  L <- length(compartments@markers)
  mult <- round(nPhi/L)
  nPhi <- L*mult
# force some phi-candidates to live near each vertex
  oneper <- max(round(nPhi/100), 10)
  point <- rep(0.02, 5)  # magic numbers: c(0.92, 0.02, 0.02, 0.02, 0.02)
  for (i in 1:5) { # do this for each compartment
    params <- point
    params[i] <- 0.92
    temp <- rdirichlet(oneper, params)
    if (exists('base')) {
      base <- rbind(base, temp)
    } else {
      base <- temp
    }
  }
# simulate a uniform set of phi-vectors
  phiset <- .reorderVectors(sampleSimplex(nPhi-5*oneper, 5))
  phiset <- rbind(base, phiset)
  likelihoods <- .computeLikelihoods(phiset, segmentdata, compartments, TRUE)
# locate the maximum likelihood for each phi-vector
  maxLikeIndex <- apply(likelihoods, 1, which.max)
  phipick <- phiset[maxLikeIndex,]
# return the results, bundled in an S4 class
  new("PrefitCloneModel",
      data=segmentdata,
      phiset=phiset,
      likelihoods=likelihoods,
      maxLikeIndex = maxLikeIndex,
      phipick = phipick,
      phiv=as.vector(phipick))
}

updatePhiVectors <- function(object, compartments) {
  if (!inherits(object, "PrefitCloneModel")) {
    object <- PrefitCloneModel(object, nPhi)
  }
  nPhi <- nrow(object@phiset)
  L <- nrow(object@data)
  multiplier <- round(nPhi/L)
# resample the phi vectors to be near the ones selected as
# optimal in the first pass
  newphiset <- matrix(NA, ncol=5, nrow=nPhi)
  for (i in 1:L) {
    index <- 1 + multiplier*(i-1)
    iset <- index:(index+multiplier-1)
    newphiset[iset,] <- rdirichlet(multiplier, 2*object@phiset[object@maxLikeIndex[i],])
  }
  newphiset <- .reorderVectors(newphiset)
# get the likelihoods for the new phis
  likelihoods <- .computeLikelihoods(newphiset, object@data, compartments, TRUE)
  maxLikeIndex <- apply(likelihoods, 1, which.max)
  phipick <- newphiset[maxLikeIndex,]
  new("PrefitCloneModel",
      data=object@data,
      phiset=newphiset,
      likelihoods=likelihoods,
      maxLikeIndex = maxLikeIndex,
      phipick = phipick,
      phiv = as.vector(phipick))
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
