### R code from vignette source 'cloneFinder.Rnw'

###################################################
### code chunk number 1: lib
###################################################
library(CloneFinder)


###################################################
### code chunk number 2: commod
###################################################
showClass("CompartmentModel")


###################################################
### code chunk number 3: compartments
###################################################
set.seed(2726642) # for reproducible examples
nSeg <-  1000     # number of segments supposedly found by CBS
markers <- round(runif(nSeg, 25, 1000))  # numbers of markers
# set 'known' centers for the pure compartments
xy <- data.frame(x = c(0.2, 0.7, 0.8, 0.1, 0.4),
                 y = c(0.2, 0.3, 0.5, 0.9, 0.7))
# build the model. sigma0 = std dev at one marker
baseModel <- CompartmentModel(markers, xy, sigma0=0.25)
rm(nSeg, xy)


###################################################
### code chunk number 4: tumor
###################################################
wts <- rev(5^(1:5))
wts <- wts/sum(wts)  # prevalence of differnt compartments
fracs <- c(5, 3, 1)  # relative frequency of subclones
# length of 'fracs' is the number of clones
# now simulate a tumor;
tumor <- Tumor(baseModel, fracs, wts)
rm(wts, fracs, markers, baseModel)
class(tumor)


###################################################
### code chunk number 5: simdata
###################################################
simdata <- generateData(tumor)


###################################################
### code chunk number 6: cloneFinder.Rnw:108-109
###################################################
sizeplot(simdata, tumor)


###################################################
### code chunk number 7: pcm
###################################################
pcm <- PrefitCloneModel(simdata)
summary(pcm)


###################################################
### code chunk number 8: cloneFinder.Rnw:170-171
###################################################
hist(pcm, breaks=77)


###################################################
### code chunk number 9: upv
###################################################
upv <- updatePhiVectors(pcm)


###################################################
### code chunk number 10: cloneFinder.Rnw:211-225
###################################################
hist(upv, breaks=55)

lp <- function(p) log(p/(1-p))
ea <- function(a) {
  temp <- exp(a)
  temp/(1+temp)
}

dang <- lp(upv@phipick)
L <- 25
dang[dang < -L] <- NA
dang[dang > L] <- NA
hist(dang, breaks=55)



