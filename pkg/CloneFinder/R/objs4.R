########################################################################
# Part 1: Simulating Data

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
Clone <- function(nSegments, weights=rep(1/5, 5)) {
  # nSegments = integer, the number of segments
  # weights   = vector, the prevalence of each compartment
  
  # start with sanity checks on the weights
  if(any(is.na(weights))) stop("Missing weights.")
  if(any(weights < 0)) stop("Negative weights.")
  weights <- weights/sum(weights) # ensure they sum to 1

  # now sample the compartments to generate a clone
  segs <- sample(length(weights), nSegments, replace=TRUE, prob=weights)
  names(segs) <- paste("Segment", 1:nSegments, sep='')
  new("Clone", segments=segs, weights=weights)
}

############ TUMOR ############
# A tumor is a set of clones, each of which is associated with a
# fraction, subject to the constraint that the sum of the fractions
# equals one.
setClass("AbstractTumor", representation=list(
                            data = "matrix",
                            fraction = "numeric",
                            markers = "numeric",
                            weights = "numeric"
                            ))

AbstractTumor <- function(fracs, markers, weights) {
  # fracs   = vector, the fraction of cells represented by each clone
  # markers = vector, the number of SNP markers per segment
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

  # we allow the user to pass in 'nSegments' instead of the vector
  # of markers. In that case, we simulate them internally
  if (length(markers) == 1) {
    markers <- round(runif(markers, 25, 2000))
  }
  nSegments <- length(markers)

  # fill in the matrix of compartments, one clone-column at a time
  temp <- matrix(NA, nrow=nSegments, ncol=L)
  for (i in 1:L) {
    temp[,i] <- s <- Clone(nSegments, weights)@segments
  }
  dimnames(temp) <- list(names(s), names(fracs))
  if(is.null(names(markers))) names(markers) <- names(s)
  new("AbstractTumor", data = temp, fraction=fracs, markers=markers, weights=weights)
}

############ TUMOR REPRESENTATION ############
setClass("Tumor", 
         contains="AbstractTumor",
         slots=c(
           compartments = "matrix",
           pureCenters = "data.frame",
           centers = "data.frame"
           )
         )

# We should be able to convert the "segments x clone" matrix representation
# of a tumor into a "segments x compartments" representation
#
Tumor <- function(object, xy) {
  # object = Tumor produced by the 'genTumor' function
  # xy     = matrix containing the fixed centers of compartments in columns 'x' and 'y'
   
  if (!inherits(object, "AbstractTumor")) stop("The 'object' must belong to the 'AbstractTumor' class.")
    nCompartments <- length(object@weights)
    nSegments <- nrow(object@data)
    fvec <- matrix(object@fraction, ncol=1)
    repr <- matrix(NA, nrow=nSegments, ncol=nCompartments)
    for (comp in 1:nCompartments) {
        temp <- 1*(object@data == comp)
        repr[,comp] <- as.vector(temp %*% fvec)
    }
    dimnames(repr) <- list(rownames(object@data), names(object@weights))
    centers <- as.data.frame(repr %*% as.matrix(xy))
    new("Tumor", object,
        compartments=repr, pureCenters=xy, centers=centers)
}

############ SIMULATING DATASET ############

# input:
#   'centers; is the result of calling "findCenters"
#   'markers' is the vector with the number of markers per segment
generateData <- function(object, sigma0=0.25) {
  # object = result of calling "representTumor"
  if (!inherits(object, "Tumor"))
    stop(paste("Incorrect class of 'object':", class(object)))
  centers <- object@centers
  markers <- object@markers
  sigma <- sigma0/sqrt(markers) # standard error of the mean
  xvec <- rnorm(length(markers), centers$x, sigma)
  yvec <- rnorm(length(markers), centers$y, sigma)
  data.frame(x=xvec, y=yvec)
}

########################################################################
# PArt 2: Bayesian Model Assessment

# compute the likelihood for each row in the dataset
#
# TODO: Understand the role of sigma0.
likely <- function(dataset, phi, object, sigma0=1) {
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
  sigma <- sigma0/sqrt(markers) # SEM
  px <- dnorm(dataset$x, center$x, sigma)
  py <- dnorm(dataset$y, center$y, sigma)
  px*py
}


sampleSimplex <- function(n, d=5) {
    result <- matrix(NA, nrow=n, ncol=d)
    for (i in 1:n) {
        result[i,] <- diff(sort( c(0, 1, runif(4, 0, 1)) ))
    }
    result
}

####################################################
# from GlobalMaxLike

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
