library(CloneFinder)

#############
# check that 'forward' works for both vectors and matrices
# and that it returns a value of the same type
mat <- matrix(c(0.6, 0.3, 0.1, 0.2, 0.3, 0.5), byrow=TRUE, ncol=3)
vec <- mat[1,]
CloneFinder:::forward(mat)
CloneFinder:::forward(vec)
class(CloneFinder:::forward(vec))
rm(mat, vec)

# check forward and backward for both vectors and matrices
# matrix input
mat <- matrix(c(0.6, 0.3, 0.1, 0.2, 0.3, 0.5), byrow=TRUE, ncol=3)
fm <- CloneFinder:::forward(mat)
bm <- CloneFinder:::backward(fm)
all(round(mat - bm, 15)==0) # round-trip agrees to 15 decimal places

# vector input
vec <- mat[1,]
fv <- CloneFinder:::forward(vec)
bv <- CloneFinder:::backward(fv)
class(bv)
all(round(vec - bv, 15)==0) # round-trip agrees to 15 decimal places

# bigger matrix
x <- runif(10000, 0, 1)
V <- matrix(CloneFinder:::lp(x), ncol=4)
F <- CloneFinder:::backward(V) # nonneg and sum to 1
temp <- apply(F, 1, sum)
summary(temp)

w <- CloneFinder:::forward(F)
flap <- CloneFinder:::backward(w)
summary(as.vector(F - flap))
rm(mat, fm, bm, vec, fv, bv, x, V, F, temp, w, flap)

#############
# check vectorization formula
Zs <- array(1:60, dim=c(4, 5, 3))
psi <- 1:3
y <- sweep(Zs, 3, psi, "*")
phinew <- apply(y, 1:2, sum)
phinew
rm(Zs, psi, y, phinew)

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
pcm <- PrefitCloneModel(dataset)
# update it
upd <- updatePhiVectors(pcm)
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
