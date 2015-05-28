########################################################################
# Part 1: Bayesian Model Assessment, phi-vectors per segment
#

###### COMPARTMENT MODEL ######
# The underlying idea is that there is a set of pure "compartments"
# representing fundamental copy number states.  In the current model,
# these are characterized by their centers.  The variability of the
# observations, however, depends on (1) the number of markers in each
# segment and (2) an estimate, sigma0, of the core variability in
# measurements observed at a single marker.

setClass("CompartmentModel",
         representation=list(
             markers = "numeric",             
             pureCenters = "data.frame",
             sigma0="numeric"))

CompartmentModel <- function(markers, pureCenters, sigma0) {
  # we allow the user to pass in 'nSegments' instead of the vector
  # of markers. In that case, we simulate them internally
  if (length(markers) == 1) {
    markers <- round(runif(markers, 25, 2000))
  }
  nSegments <- length(markers)
  if(is.null(names(markers))) names(markers) <- paste("Segment", 1:nSegments, sep='')

  new("CompartmentModel",
    markers=markers,
    pureCenters=pureCenters,
    sigma0=sigma0)
}

###### Likelihoods computations #####
# Given the compartment model (which contains the centers
# and the number of markers per segment) and the vector of
# fraction per compartment, 
# compute the likelihood for each row in the dataset
likely <- function(dataset, phi, compartmentModel, log=FALSE) {
  # dataset = matrix with 'x' and 'y' columns produced by "generateData"
  # phi = vector of probabilities for each compartment
  xy <- compartmentModel@pureCenters
  markers <- compartmentModel@markers
  sigma0 <- compartmentModel@sigma0
  if (length(phi) < nrow(xy)-1) stop("You did not supply enough entries in 'phi'")
  if (length(phi) > nrow(xy)) stop("You supplied too many entries in 'phi'")
  if (length(phi) < nrow(xy)) {
    lastphi <- 1 - sum(phi)
    phi <- c(phi, lastphi)
  }
  if (any(phi < 0)) stop("Negative probabilities are not allowed")
  phi <- matrix(phi/sum(phi), nrow=1) # make sure they add up to 1
  center <- as.data.frame(phi %*% as.matrix(xy))
  secondMoment <- phi %*% (xy^2 + sigma0^2)
  sigma <- sqrt(sweep(secondMoment - center^2, 1, markers, '/'))
  px <- dnorm(dataset$x, center$x, sigma[,1], log)
  py <- dnorm(dataset$y, center$y, sigma[,2], log)
  if(log) {
      value <- px + py
  } else {
      value <- px * py
  }
  value
}

sampleSimplex <- function(n, d=5) {
    result <- matrix(NA, nrow=n, ncol=d)
    for (i in 1:n) {
        result[i,] <- diff(sort( c(0, 1, runif(d-1, 0, 1)) ))
    }
    result
}

########################################################################
# Part 2: Simulating Data

############ CLONE ############
# A clone is a list of segments, where each segment represents exactly
# one pure compartment.

setClass("Clone", representation=list(
                    segments = "integer",
                    weights = "numeric"
                    ))

# We should be able to simulate a clone if we know the number of segments,
# the number of compartments, and the relative frequency/prevalence of
# each compartment.
Clone <- function(nSegments, weights=rep(1/5, 5), segnames=NULL) {
  # nSegments = integer, the number of segments
  # weights   = vector, the prevalence of each compartment

  if (nSegments < 1) stop("Number of segments must be positive.")
  
  # start with sanity checks on the weights
  if(any(is.na(weights))) stop("Missing weights.")
  if(any(weights < 0)) stop("Negative weights.")
  weights <- weights/sum(weights) # ensure they sum to 1

  # now sample the compartments to generate a clone
  segs <- sample(length(weights), nSegments, replace=TRUE, prob=weights)
  if (is.null(segnames)) segnames <- paste("Segment", 1:nSegments, sep='')
  if (length(segnames) != length(segs)) stop("Wrong number of segment names.")
  names(segs) <- segnames
  new("Clone", segments=segs, weights=weights)
}

############ TUMOR ############
# A tumor is a set of clones, each of which is associated with a
# fraction, subject to the constraint that the sum of the fractions
# equals one.
setClass("Tumor",
         contains = "CompartmentModel",
         slots=list(
             data = "matrix",         # segments x clones
             fraction = "numeric",
             weights = "numeric",
             compartments = "matrix", # segments x compartmentys
             centers = "data.frame",
             SEM="matrix"
             ))

Tumor <- function(object, fracs, weights) {
  # object = CompartmentModel
  # fracs   = vector, the fraction of cells represented by each clone
  # weights = vector, the prevalence of each pure compartment

  # sanity checks on the fractions
  if(any(is.na(fracs))) stop("Missing clone fraction.")
  if(any(fracs < 0)) stop("Negative clone fraction.")
  fracs <- fracs/sum(fracs) # ensure they sum to 1
  L <- length(fracs)
  if(is.null(names(fracs))) names(fracs) <- paste("Clone", 1:L, sep="")

  # we could do the same check on the weights, but they are just
  # being passed on to 'Clone' which already checks them
  if (is.null(names(weights)))
    names(weights) <- paste("Compartment", 1:length(weights), sep="")

  nSegments <- length(object@markers)
  segnames <- names(object@markers)

  # fill in the matrix of compartments, one clone-column at a time
  dataset <- matrix(NA, nrow=nSegments, ncol=L)
  for (i in 1:L) {
    dataset[,i] <- s <- Clone(nSegments, weights, segnames)@segments
  }
  dimnames(dataset) <- list(names(s), names(fracs))

  # convert from the "segments x clones" representations in 'dataset'
  # to a "segments x compartments" representation
  nCompartments <- length(weights)
  fvec <- matrix(fracs, ncol=1)
  repr <- matrix(NA, nrow=nSegments, ncol=nCompartments)
  for (comp in 1:nCompartments) {
      temp <- 1*(dataset == comp)
      repr[,comp] <- as.vector(temp %*% fvec)
  }
  dimnames(repr) <- list(rownames(dataset), names(weights))
  # get the averaged centers, over compartments
  xy <- as.matrix(object@pureCenters)
  centers <- as.data.frame(repr %*% xy)
  # get the SEM
  secondMoment <- repr %*% (xy^2 + object@sigma0^2)
  sigma <- sqrt(sweep(secondMoment - centers^2, 1, object@markers, '/'))
  
  new("Tumor", object,
      data = dataset,
      fraction=fracs,
      weights=weights,
      compartments=repr,
      centers=centers,
      SEM=sigma
      )
}

trueZ <- function(tumor) {
  Zmat <- array(0, dim=c(nrow(tumor@centers), # segments
                     nrow(tumor@pureCenters), # compartments
                     ncol(tumor@data)))       # clones
  for (wclone in 1:ncol(tumor@data)) {
    for(wseg in 1:nrow(tumor@centers)) {
      wcomp <- tumor@data[wseg, wclone]
      Zmat[wseg,wcomp,wclone] <- 1
    }
  }
  Zmat
}

setMethod("summary", signature("Tumor"), function(object, ...) {
    tdata <- object@data
    table(A=tdata[,1], B=tdata[,2])
})

############ SIMULATING DATASET ############

# input:
#   'centers; is the result of calling "findCenters"
#   'markers' is the vector with the number of markers per segment
generateData <- function(object) {
  if (!inherits(object, "Tumor"))
    stop(paste("Incorrect class of 'object':", class(object)))
  centers <- object@centers
  markers <- object@markers
  xy <- object@pureCenters
  # now we can genrate the data
  xvec <- rnorm(length(markers), centers$x, object@SEM[,1])
  yvec <- rnorm(length(markers), centers$y, object@SEM[,2])
  data.frame(x=xvec, y=yvec)
}

sizeplot <- function(simdata, tumor) {
    size <- 1 + round(sqrt(tumor@markers)/15)/2
    filt <- size > log(128)/3
    x <- tumor@pureCenters$x
    y <- tumor@pureCenters$y
    eps <- 0.1*diff(range(x))
    plot(x, y, type='n',
         xlim=c(min(x)-eps, max(x)+eps),
         ylim=c(min(y)-eps, max(y)+eps))
    text(x, y, 1:nrow(tumor@pureCenters))
    points(tumor@pureCenters, cex=15, col='gray25')
    points(simdata$x[filt], simdata$y[filt],
           cex=size[filt], col="gray35", pch=16)
    invisible(simdata)
}


####################################################
# from GlobalMaxLike
#  KRC: DO we ened these any more? Or are they part of an
# old approach that didn't generalize?

makeLatent <- function(vector){
  clone <- matrix(0, ncol=5, nrow=length(vector))
  for (i in 1:length(vector)){
    clone[i, vector[i]] <- 1
  }
  clone
}

estBetaParams <- function(vec) {
  mu <- mean(vec)
  s2 <- var(vec)
  temp <- mu*(1-mu)/s2 - 1
  alpha <- mu*temp
  beta <- (1 - mu)*temp
  c(alpha = alpha, beta = beta)
}


maxfun <- function(segment){
  d <- temp[segment,]
  range <- 1:1000
  range <- range[!range == segment]
  optimfun <- function(n){
    distfun <- function(segment){
      log(dbeta(d[n], paramtable[segment, 1], paramtable[segment, 2]))
    }
    sum(sapply(range, distfun))+log(d[n]/sum(d))
  }
  d <- sapply(1:10000, optimfun)
  i <- which.max(d)
  newphiset[i,]
}
