source("objs4.R")

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

# example: tumor with two clones of unequal weight
tumor <- AbstractTumor(c(3,1), 200, wts)
tdata <- tumor@data
summary(tdata)
table(A=tdata[,1], B=tdata[,2])
summary(tumor@markers)
tumor@fraction
rm(tumor, tdata)

# example: specifying segment lengths
segs <- runif(300, 20, 3000)
tumor <- AbstractTumor(c(3,1), segs, wts)
rm(segs)
tdata <- tumor@data
table(A=tdata[,1], B=tdata[,2])
summary(tumor@markers)
tumor@fraction
tumor@weights

# continued example
# known compartment centers
xy <- data.frame(x = c(.2, .7, .8, .1, .4),
                 y = c(.2, .3, .5, .9, .7))
temp <- Tumor(tumor, xy)
head(tumor@data, 10)
head(temp@compartments, 10)
rm(temp, tumor, tdata, xy, wts)

