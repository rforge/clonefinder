e4 <- guessPsi(upd, 4)
e4
# refine with EM-algorithm
f4 <- runEMalg(e4, dataset, tumor)
f4$psi
f4$loglike

e5 <- guessPsi(upd, 5)
e5
# refine with EM-algorithm
f5 <- runEMalg(e5, dataset, tumor)
f5$psi
f5$loglike

e2 <- guessPsi(upd, 2)
e2
# refine with EM-algorithm
f2 <- runEMalg(e2, dataset, tumor)
f2$psi
f2$loglike

e1 <- 1
f1 <- runEMalg(e1, dataset, tumor)

park <- 0*upd@phipick
for (i in 1:nrow(park)) {
  w <- which.max(upd@phipick[i,])
  park[i,w] <- 1
}

loglikes <- sum(tock <- sapply(1:nrow(park), function(i, phi) {
    ci <- compModel
    ci@markers <- ci@markers[i]
    sum(log(likely(dataset[i,,drop=FALSE], phi[i,], ci)))
  }, phi=park))
loglikes # really ?!?

ll <- c(f2$loglike, final$loglike, f4$loglike, f5$loglike)

# assume pr(n) is the prior that the number of true clones is n
# then the posterior po(n) = pr(n)*like(data|n) / like(data)
# so log(po(n)) = log(pr(n)) + loglike(data|N) - loglike(data)

