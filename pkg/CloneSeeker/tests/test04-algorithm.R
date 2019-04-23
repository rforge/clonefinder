if (packageVersion("CloneSeeker") < "0.9.0") {
  stop("You need to update 'CloneSeeker'.")
}
library("CloneSeeker")

psiSets <- list(c(1),                     # only one clone
                c(0.7, 0.3),              # two clones
                c(0.6, 0.25, 0.15),       # three clones
                c(0.5, 0.25, 0.15, 0.10)) # four clones

paramSets <- data.frame(nu = rep(c(0, 100, 100), each = 2),
                        pcnv = rep(c(1, 0.5, 0), each = 2),
                        norm.contam = rep(c(TRUE, FALSE), times = 3))

paramSets <- data.frame(nu = rep(c(0, 100, 100)),
                        pcnv = rep(c(1, 0.5, 0)))

#########################################
### Simulate the data sets.
RNGversion("3.5.3")
set.seed(20248)
simData <- list()
for (psi in psiSets) {
  cat("Working on psi =", psi, "\n", file=stderr())
  for (J in 1:nrow(paramSets)) {
    cat("Working on set J =", J, "\n", file=stderr())
    tumor <- Tumor(psi,
                   rounds = 100,
                   cnmax = 4,
                   nu = paramSets$nu[J],
                   pcnv = paramSets$pcnv[J],
                   norm.contam = FALSE)
    dset <- generateTumorData(tumor,
                              snps.seq = 10000,
                              snps.cgh = 600000,
                              mu = 70,
                              sigma.reads = 25,
                              sigma0.lrr = 0.15, 
                              sigma0.baf = 0.03,
                              density.sigma = 0.1)
    simData[[1 + length(simData)]] <- list(tumor = tumor,
                                           dset = dset,
                                           psi = psi,
                                           params = paramSets[J,])
  }
}
length(simData)

#########################################
### Set up algorithm parameters
psis.20 <- generateSimplex(20,5)
cnmodels <- as.matrix(expand.grid(lapply(1:5, function(i){ 0:5 })))
pars <- list(sigma0=5, theta = 0.9, ktheta = 0.3, mtheta = 0.9,
             alpha = 0.5, thresh = 0.04, cutoff = 100, Q = 100, iters = 4)

### Don't run full tests at CRAN.
if (Sys.getenv("KRC_DO_ALL_TESTS") == "TRUE") {
  testset <- 1:length(simData)
} else {
  testset <- c(1, 5)
}

for (J in testset) {
  cat("\n\nDataset", J, "\n", file=stdout())
  dset <- simData[[J]]$dset
  ## Cheat and limit the number of mutations to speed up the algorithm for CRAN
  ## This code was copied out of the "seekClones" algorithm.
  seqdata <- dset$seq.data
  if(nrow(seqdata) > 0) {
    read.den <- density(seqdata$refCounts)
    peak <- read.den$x[which.max(read.den$y)]
  }
  mutdata <- seqdata[seqdata$status == 'somatic',]
  if (nrow(mutdata) > 1) {
    mut.filt <- CloneSeeker:::filterMutations(mutdata, mu=peak, threshold=3)
    seqdata <- mut.filt$mat
    if (nrow(seqdata) > 10) seqdata <- seqdata[1:10,]
  }
  ## "Try" was placed here during debugging. These should succeed.
  ra <- try( runAlg(dset$cn.data, seqdata,
                    cnmodels, psis.20,
                    pars = pars, imputedCN = NULL) )
  if (inherits(ra, "try-error")) {
    cat(ra, "\n", stdout())
  } else {
    cat("psi:\n")
    print(ra$psi)
    cat("A,B:\n")
    print(cbind(A = ra$A, B = ra$B))
    cat("eta:\n")
    print(data.frame(etaA = ra$etaA, etaB = ra$etaB))
    cat("indices:\n")
    print(lapply(ra$indices, summary))
    cat("filtered data:\n")
    print(lapply(ra$filtered.data, summary))
    print(summary(ra$mutated))
    cat("posteriors\n")
    print(summary(ra$psiPosts))
  }
}

if (FALSE) { # save the profiling code for later
  Rprof("profile3.txt")
  for (J in c(1, 5, 9)) {
    dset <- simData[[J]]$dset
    ra <- seekClones(dset$cn.data, dset$seq.data,
                     cnmodels, psis.20,
                     pars = pars, imputedCN = NULL)
  }
  Rprof(NULL)
  summaryRprof("profile3.txt")
}
