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

PrefitCloneModel <- function(segmentdata, compartments, nPhi = 20000) {
# use the Jeffreys (Dirichlet) prior on phi
  N <- nrow(compartments@pureCenters)
  phiset <- rdirichlet(nPhi, rep(1/2, N))
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

# posterior distributions on phi
posteriorPhi <- function(pcm) {
  t(apply(pcm@likelihoods, 1, function(x) {
    maxx <- max(x)
    total <- sum(ptwise <- exp(x - maxx))
    ptwise/total
  }))
}
#  euclidean distance in phi-space
ed <- function(pcm, phi0) {
  sqrt(apply(sweep(pcm@phiset, 2, phi0, "-")^2, 1, sum))
}

plotPosteriorCDF <- function(pcm, phi0, irow, post=NULL, ...) {
  if(is.null(post)) post <- posteriorPhi(pcm)
  delta <- ed(pcm, phi0)
  od <- order(delta)
  cdf <- cumsum(post[irow, od]) # cumulative distribution by radius
  plot(sort(delta), cdf, xlab="Radius", ylab="CDF", ...)
  invisible(pcm)
}

posteriorRegion <- function(pcm, radius=0.05, phi0=NULL, post=NULL) {
  if(is.null(post)) post <- posteriorPhi(pcm)
  targets <- t(sapply(1:nrow(post), function(i, phi0, rset) {
    if(is.null(phi0)) phi0 <- pcm@phiset[pcm@maxLikeIndex[i],]
    delta <- ed(pcm, phi0)
    od <- order(delta)
    sapply(rset, function(r) {
      n <- which(delta[od] >= r)[1] - 1
      sum(post[i,od[1:n]])
    })
  }, phi0=phi0, rset=radius))
  colnames(targets) <- paste("R", radius, sep='')
  targets
}

posteriorQuantile <- function(pcm, q=c(0.7, 0.8, 0.9, 0.95, 0.99), phi0=NULL, post=NULL) {
  if(is.null(post)) post <- posteriorPhi(pcm)
  targets <- t(sapply(1:nrow(post), function(i, phi0, qset) {
    if(is.null(phi0)) phi0 <- pcm@phiset[pcm@maxLikeIndex[i],]
    delta <- ed(pcm, phi0)
    od <- order(delta)
    cs <- cumsum(post[i,od])
    sapply(qset, function(q) {
      delta[od][which(cs>q)][1]
    })
  }, phi0=phi0, qset=q))
  colnames(targets) <- paste("Q", q, sep='')
  targets
}
