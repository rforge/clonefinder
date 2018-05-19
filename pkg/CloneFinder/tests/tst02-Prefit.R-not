library(CloneFinder)

set.seed(363453)

# generate the markers explicitly
nSeg <- 1000
markers <- round(runif(nSeg, 25, 1000))
# Centers needed to define the pure structure
xy <- data.frame(x = c(.2, .7, .8, .1, .4),
                 y = c(.2, .3, .5, .9, .7))
# start by creating the compartment model
baseModel <- CompartmentModel(markers, xy, sigma0=0.25)
rm(xy, nSeg, markers)

# Now we set up an abstract tumnor
wts <- rev(5^(1:5))
wts <- wts/sum(wts)
fracs <- c(5, 3, 1)
TrueNclones <- length(fracs)
# now simulate a tumor; length of 'fracs' in first argument is number of clones
tumor <- Tumor(baseModel, fracs, wts)
rm(wts, fracs)
ls()

# simulate data by selecting the weighted means with appropriate
# standard error of the mean
simdata <- generateData(tumor)

firstPass <- PrefitCloneModel(simdata, baseModel)
plot(firstPass)
hist(firstPass, breaks=123)
summary(firstPass)
