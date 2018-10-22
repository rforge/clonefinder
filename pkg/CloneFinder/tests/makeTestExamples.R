if (packageVersion("CloneFinder") < "0.8.7") {
  stop("You need to update 'CloneFinder'.")
}
library("CloneFinder")

psiSets <- list(c(1),                     # only one clone
                c(0.7, 0.3),              # two clones
                c(0.6, 0.25, 0.15),       # three clones
                c(0.5, 0.25, 0.15, 0.10)) # four clones

paramSets <- data.frame(nu = rep(c(0, 100, 100), each = 2),
                        pcnv = rep(c(1, 0.5, 0), each = 2),
                        norm.contam = rep(c(TRUE, FALSE), times = 3))

paramSets <- data.frame(nu = rep(c(0, 100, 100)),
                        pcnv = rep(c(1, 0.5, 0)))

set.seed(20248)
simData <- list()
for (psi in psiSets) {
  cat("Working on psi =", psi, "\n", file=stderr())
  for (J in 1:nrow(paramSets)) {
    cat("Working on set J =", J, "\n", file=stderr())
    tumor <- Tumor(psi, rounds = 400,
                   cnmax = 4,
                   nu = paramSets$nu[J],
                   pcnv = paramSets$pcnv[J],
                   norm.contam = FALSE)
    dset <- generateTumorData(tumor,
                              snps.seq = 1000000,
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
### Restart here
psis.20 <- generateSimplex(20,5)
cnmodels <- as.matrix(expand.grid(lapply(1:5, function(i){ 0:5 })))
pars <- list(sigma0=5, theta = 0.9, ktheta = 0.3, mtheta = 0.9,
             alpha = 0.5, thresh = 0.04, cutoff = 100, Q = 100, iters = 4)

K <- 5

for (J in 1:length(simData)) {
  cat("Dataset", J, "\n", file=stdout())
  dset <- simData[[J]]$dset
  ra <- try( runAlg(dset$cn.data, dset$seq.data,
                    cnmodels, psis.20,
                    pars = pars, imputedCN = NULL) )
  if (inherits(ra, "try-error")) {
    cat(ra, "\n", stdout())
  } else {
    print(ra$psi)
    print(ra$indices)
    print(ra$filtered.data)
    print(ra$etaA)
    print(ra$etaB)
    print(ra$A)
    print(ra$B)
    print(ra$mutated)
    print(summary(ra$psiPosts))
  }
}



if (FALSE) {
for (psi in psiSets) {
  cat("Working on psi =", psi, "\n", file=stdout())
  for (J in 1:nrow(paramSets)) {
    print(paramSets[J,])
    tumor <- Tumor(psi, rounds = 400,
                   cnmax = 4,
                   nu = paramSets$nu[J],
                   pcnv = paramSets$pcnv[J],
                   norm.contam = paramSets$norm.contam[J])
    seq.clones <- lapply(1:length(tumor@clones), function(i){tumor@clones[[i]]$seq})
    print(out <- class(seq.clones[[1]]))
    temp <- data.frame(paramSets[J,], nClones = length(psi), Class=out)
    if (exists("results")) {
      results <- rbind(results, temp)
    } else {
      results <- temp
    }
  }
}
results[order(results$Class),]
}
