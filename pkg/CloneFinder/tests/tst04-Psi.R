library(CloneFinder)

#############################
# simulate a sample data set
set.seed(539121)
# pure centers
xy <- data.frame(x = log10(c(2, 2, 1, 3, 4)/2),
                 y = c(1/2, 0, 0, 1/3, 1/4))
plot(xy, pch=16)
# number of segments
nSeg <- 1000
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
# update it
upd <- updatePhiVectors(pcm, tumor)
# good guess at psi-vector
estpsi <- guessPsi(upd, 3) # 3 is number of clones we aqre trying to fit
estpsi
# refine with EM-algorithm
final <- runEMalg(estpsi, dataset, tumor)
final$psi
final$loglike

Zed <- trueZ(tumor)
align <- Zed - final$Zmats

gotit <- apply(align, 1, function(x) all(x==0))
mean(gotit)
