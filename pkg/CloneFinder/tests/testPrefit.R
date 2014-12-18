library(CloneFinder)

set.seed(363453)

# Parameters that we need to define the structure
TrueNclones<- 3
nSeg <- 1000
wts <- rev(5^(1:5))
wts <- wts/sum(wts)
xy <- data.frame(x = c(.2, .7, .8, .1, .4),
                 y = c(.2, .3, .5, .9, .7))

# generate the markers explicitly
markers <- round(runif(nSeg, 25, 1000))
fracs<-c(5, 3, 1)
TrueNclones<- length(fracs)
#nclones is the number of clones we beleive there to be
nclones<- 3

# now simulate a tumor; length of 'fracs' in first argument is number of clones
abstractTumor <- AbstractTumor(fracs, markers, wts)
# and get the concrete representation
tumor <- Tumor(abstractTumor, xy)
# clean up by removing stuff we don't need
rm(markers, wts, abstractTumor)
ls()

# simulate data by selecting the weighted means with appropriate
# standard error of the mean
simdata <- generateData(tumor)

firstPass <- PrefitCloneModel(simdata)
plot(firstPass)
hist(firstPass, breaks=123)
summary(firstPass)

secondPass <- updatePhiVectors(firstPass)
plot(secondPass)
hist(secondPass, breaks=56)
summary(secondPass)
