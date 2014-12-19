library(CloneFinder)
set.seed(523476)

# example; equal prevalences of compartments
clone <- Clone(10000)
clone@weights
table(clone@segments)
rm(clone)

# example; highly biased prevalences of compartments
wts <- rev(5^(1:5))
wts <- wts/sum(wts)
clone <- Clone(10000, wts)
clone@weights
table(clone@segments)
rm(clone)

# known compartment centers
xy <- data.frame(x = c(.2, .7, .8, .1, .4),
                 y = c(.2, .3, .5, .9, .7))
baseModel <- CompartmentModel(200, xy, 0.25)
# example: tumor with two clones of unequal weight
tumor <- Tumor(baseModel, c(3,1), wts)
tdata <- tumor@data
summary(tdata)
table(A=tdata[,1], B=tdata[,2])
summary(tumor@markers)
tumor@fraction
rm(tumor, tdata)

# example: specifying segment lengths
segs <- runif(300, 20, 3000)
baseModel <- CompartmentModel(segs, xy, 0.25)
tumor <- Tumor(baseModel, c(3,1), wts)
rm(segs)
tdata <- tumor@data
table(A=tdata[,1], B=tdata[,2])
summary(tumor@markers)
tumor@fraction
tumor@weights
head(tumor@data, 10)
head(tumor@compartments, 10)
rm(tumor, tdata, xy, wts)

