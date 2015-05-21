# Computes the NEGATIVE log-likelihood for some value of psi,
# assuming that we already have the z1 and z2 matrices.
# Note that the input is assumed to be given by
#      x = lp(psi)
# and Zs = 3D array of z-matrices, so should have dimension
#     nSegments x nCompartments x nClones
#     where nClones equals length(x)
myTarget <- function(x, Zs, data, tumor) {
  if (length(x) != dim(Zs)[3]) stop('mismatched argmuent sizes')
  psi <- ea(x)
  # next formula is phi = sum_i (psi_i * Z_i)
  phinew <- apply(sweep(Zs, 3, psi, "*"), 1:2, sum)
  loglikes <- sum(tock <- sapply(1:nrow(phinew), function(i, phi) {
    sum(log(likely(data[i,], phi[i,], tumor, sigma0=0.25)))
  }, phi=phinew))
  - loglikes  # negate
}

# check algorithm
# Zs <- array(1:60, dim=c(4, 5, 3))
# psi <- 1:3
# y <- sweep(Zs, 3, psi, "*")
# phinew <- apply(y, 1:2, sum)

# TODO: [1] Replace the call to "optimize" by one to "optim".
# TODO: [2] Replace the nested loop over i-j by a single loop
# over all possible Z-assignments.  Theoretically, you could precompute
# all possible assignments (at most 5^5 = 3125) and store them in a
# matrix or array.

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
