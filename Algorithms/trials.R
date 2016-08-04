# potential replacement for PrefitCloneModel 
PCM <- function(segmentdata, compartments, nPhi = 20000) {
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

# new function to compute posterior distribution
posterior <- function(pcm) {
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

posteriorRegion <- function(pcm, radius=0.05, phi0=NULL) {
  post <- posterior(pcm)
  targets <- t(sapply(1:nrow(post), function(i, phi0, qset) {
    if(is.null(phi0)) phi0 <- pcm@phiset[pcm@maxLikeIndex[i],]
    delta <- ed(pcm, phi0)
    od <- order(delta)
    n <- which(delta[od] >= radius)[1] - 1
    sum(post[i,od[1:n]])
  }, phi0=phi0, qset=q))
  as.vector(targets)
}

posteriorQuantile <- function(pcm, q=c(0.7, 0.8, 0.9, 0.95, 0.99), phi0=NULL) {
  post <- posterior(pcm)
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

