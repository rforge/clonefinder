source("../tests/tst04-Psi.R")

e1 <- 1
f1 <- runEMalg(e1, dataset, tumor)
f1$psi
f1$loglike

e2 <- guessPsi(upd, 2)
e2
# refine with EM-algorithm
f2 <- runEMalg(e2, dataset, tumor)
f2$psi
f2$loglike

e3 <- estpsi
f3 <- final

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

e6 <- guessPsi(upd, 6)
e6
# refine with EM-algorithm
f6 <- runEMalg(e6, dataset, tumor)
f6$psi
f6$loglike


ll <- c(f1$loglike, f2$loglike, f3$loglike, f4$loglike, f5$loglike, f6$loglike)
plot(ll)

# assume pr(n) is the prior that the number of true clones is n
# then the posterior po(n) = pr(n)*like(data|n) / like(data)
# so log(po(n)) = log(pr(n)) + loglike(data|N) - loglike(data)

# We want to take an exponentially decaying prior of the form
#    Prior(d) = exp(-theta*c*d)
# The trickiest part is figuring out how to choose c and how to
# interpret this stuff sensibly
#
# Note that, by the usual rules of exponents,
#    Prior(d) = exp(-theta*d)^c
# so we can think of 'c' as a power on some simpler prior
# Assuming that each segment is independent, putting such a prior
# on each segment, one would get this overall weighted prior by
# taking c = number of segments, which for our example is 1000.

simplePrior <- function(theta, nseg) {
  d <- 1:6
  exp(-nseg*d*theta)
}
barplot(simplePrior(1, 1))

# Now we use the Auer-Gervini trick to obtain a step function
# for maximum a posteriori number of clones as a function of
# the decay parametr theta

sf <- stepfun( rev(diff(ll)) / 1000, 6:1)
plot(sf)

knots(sf)
# So, when theta is about 0.3 and c=nseg = 1000, we should
# switch from thinking there are 4 clones down to 3 clones

# This corresponds to the following prior on a single segment
barplot(simplePrior(knots(sf)[3], 1))

# We don't swithc from 3 to 2 until about 2.15, which
# corresponds to this per-segment prior
barplot(simplePrior(knots(sf)[4], 1))

# And we don't switch from 2 to 1 until aboutr 42, whihc looks like this
barplot(simplePrior(knots(sf)[4], 5))

# We can also look at the corresponding log posteriors

logpost <- function(theta, nseg) {
  d <- 1:6
  -d*nseg*theta + ll[d]
}
barplot( logpost(knots(sf)[1], 1000), main="Log Posterior (6--5)")
barplot( logpost(knots(sf)[2], 1000), main="Log Posterior (5--4)")
barplot( logpost(knots(sf)[3], 1000), main="Log Posterior (4--3)")
barplot( logpost(knots(sf)[4], 1000), main="Log Posterior (3--2)")
barplot( logpost(knots(sf)[5], 1000), main="Log Posterior (2--1)")
