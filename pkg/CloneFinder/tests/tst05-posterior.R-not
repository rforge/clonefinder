library(CloneFinder)

#############################
# simulate a sample data set
set.seed(539121)
# pure centers
xy <- data.frame(x = log10(c(2, 2, 1, 3, 4)/2),
                 y = c(1/2, 0, 0, 1/3, 1/4))
# small number of segments
nSeg <- 10
# number of SNP markers per segment
markers <- round(runif(nSeg, 25, 1000))
compModel <- CompartmentModel(markers, xy, 0.25)
# probability of a pure clonal segment in each compartment
wts <- rev(5^(1:5))
wts <- wts/sum(wts)
# percentage of cells in each (of three) clone(s)
psis <- c(0.6, 0.3, 0.1)

tumor <- Tumor(compModel, psis, wts)
dataset <- generateData(tumor)

# prefit the model
pcm <- PrefitCloneModel(dataset, tumor)
#############################
# now test the posterior distributions

post <- posteriorPhi(pcm)
dim(post) # segments by phi-vectors
summary(apply(post, 1, sum)) # all 1's

dataset # pure normal for 1, 2, 3
# should make this an exported function
phi0 <- c(1,0,0,0,0) # euclidean distance from normal vertex
plotPosteriorCDF(pcm, phi0, irow=1, post=post, type='l')
plotPosteriorCDF(pcm, phi0, irow=8, post=post, type='l')


# radii needed to acheive a given quantile of posterior probability
postq <- posteriorQuantile(pcm, phi0=phi0, post=post)
postq

# probability of true phi being within a given radius of a phi-vector
pr <- posteriorRegion(pcm, radius=c(0.05, 0.10), phi0=phi0, post=post)
round(pr, 3)

