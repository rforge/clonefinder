source("C:\\Users\\Mark\\Documents\\Lab\\Code\\objects.R")
source("C:\\Users\\Mark\\Documents\\Lab\\Code\\PosteriorFunctionMethod2 (2).R")
v<-c(foo[[2]][,1])
w<-c(foo[[2]][,2])
length(v)
segments<-1000

matrixfun<-function(vector){
  c<- c(rep(0, length(vector)))
  clone<-data.frame(c1=c, c2=c, c3=c, c4=c, c5=c)
  for (i in 1:length(vector)){
    clone[i, vector[i]]<-1
  }
  clone
}
Z1<-matrixfun(v)
Z2<-matrixfun(w)

PhiMatrix<-foo$fraction[1]*Z1+foo$fraction[2]*Z2
PhiMatrix
#True Z matrix
#-----------------------------------------------

d2<- lapply(1:nPhi, postfunc)
class(d2)       # this is a list
length(d2)      # one entry for each phi-vector
class(d2[[1]])  # each entry is a numeric vector
length(d2[[1]]) # with one entry for each segment

temp <- t(as.data.frame(d2))
rownames(temp) <- NULL
temp[1:5, 1:5]
dim(temp)

temp <- apply(phiset, 1, function(C) {
  likely(d, mu1, mu2, mu3, mu4, mu5, C[1], C[2], C[3], C[4])
})

maxLikeIndex <- apply(temp, 1, which.max)
plot(maxLikeIndex)

phipick <- phiset[maxLikeIndex,]

#pass2
library(mc2d)
multiplier <- nPhi/segments
newphiset <- matrix(NA, ncol=5, nrow=nPhi)
for (i in 1:segments) {
  index <- 1 + multiplier*(i-1)
  iset <- index:(index+multiplier-1)
  newphiset[iset,] <- rdirichlet(multiplier, 2*phiset[maxLikeIndex[i],])
}
temp <- apply(newphiset, 1, function(C) {
  likely(d, mu1, mu2, mu3, mu4, mu5, C[1], C[2], C[3], C[4])
})

maxLikeIndex <- apply(temp, 1, which.max)
plot(maxLikeIndex)

phipick <- newphiset[maxLikeIndex,]
hist(phipick, breaks=123)
plot(density(phipick))
#-----------------------------------
#Now, computing the Z matrices from the data:
#first, we have to (?) estimate the number of clones. Here clearly two, but more rigorously?
library(TTR)
library(forecast)
d<-ma(phipick, 3, TRUE)
plot(density(na.omit(d)))
density(phipick)

#Assume we can get a rigorous estimate of the number of clones, here of course 2:
nclones<-2
head(phipick)
f<-function(row){
  d<-phipick[row,]
  y<-c(sort(d, decreasing=TRUE)[1:2])
  x<-c(which(d %in% y[1]), which(d %in% y[2]))
  c(x[1], y[1], x[2], y[2])
}
d<-t(data.frame(sapply(1:1000, f)))
head(d)
plot(d[,2], xlim=c(0, 1000), ylim=c(0, 1))
points(d[,4])
plot(density(d[,4]), xlim=c(0, 1), ylim=c(0, 10))
lines(density(d[,2]))
#Now we need to identify minima to know where to allocate compartments; at each segment, there 
#should be two possibilities: .25 and .75); and 1 and (when both clones are in te same 
#compartment). Here we can of course estimate minima of .12, .5, and .86.

data1<-d[,1:2]
data2<-d[,3:4]

g<-function(row){
if(data1[row, 2]>.87){
  clone1<-data1[row, 1]
}else if(data1[row, 2]<.13){
  clone1<-data2[row, 1]
}else{
  clone1<-data1[row, 1]
}
if(data1[row, 2]>.87){
  clone2<-data1[row, 1]
}else if(data1[row, 2]<.13){
  clone2<-data2[row, 1]
}else{
  clone2<-data2[row, 1]
}
c(clone1, clone2)
}
#the clone that consisting .75 of all cells is defined as clone 1
clones<-t(sapply(1:1000, g))
head(clones)

#One way to accurately estimate clone number may be to compare the success rate of predictions and
#see which clone estimate has the highest success rate, and see if this seems to invariably correspond
#some measure, maybe number of local maxima between the ones by 0  and 1? Need more rigorous method.

vectorfun<-function(vector){
  clone<-matrix(
    c(rep(0, 5000)),
    nrow=1000,
    ncol=5)
  for (i in 1:1000){
    clone[i, vector[i]]<-1
  }
  clone
}
clone1<-vectorfun(c(clones[,1]))
z<-clone1==Z1
length(which(z==FALSE))/2
#success rate: 99.4%

clone2<-vectorfun(c(clones[,2]))
y<-clone2==Z2
length(which(y==FALSE))/2
#success rate of 97.5%

#Let's look at the false assignments:
diff1<-foo$data[,1]==clones[,1]
diff2<-foo$data[,2]==clones[,2]
F1<-which(diff1 %in% FALSE)
F2<-which(diff2 %in% FALSE)

#Some ways of viewing the results...
h<-function(n, clone){
clones[n, clone]
}

df.false1<-t(rbind(mapply(h, F1, 1)))
df.false2<-t(rbind(mapply(h, F2, 2)))
             
i<-function(n, clone){
  foo$data[n, clone]
}

true1<-mapply(i, F1, 1)
estimate1<-df.false1[,1]
df1<-data.frame(true1, estimate1)

true2<-mapply(i, F2, 2)
estimate2<-df.false2[,1]
df2<-data.frame(true2, estimate2)

j<-function(n){
phipick[n,]
}

df.false1.phis<-t(rbind(sapply(F1, j)))
df.false2.phis<-t(rbind(sapply(F2, j)))

#the following data frames show phi values and compartment assignments for false estimates of each
# of the two clones seperately.
cbind(df1, as.data.frame(df.false1.phis))
cbind(df2, as.data.frame(df.false2.phis))
#All false assignments appear to be due to error in the data, not the algorithm.

#Now, using kmeans to assign sements to compartments:

data1.phi<-data1[,2]
clusters1<-as.vector(kmeans(data1.phi, 2, 20)[[1]])

data2.phi<-data2[,2]
clusters2<-as.vector(kmeans(data2.phi, 2, 20)[[1]])

k<-function(n){
if(clusters1[n]==2){
  comp<-c(data1[n, 1], data1[n, 1])
}
else{
  comp<-c(data1[n, 1], data2[n, 1])
}
comp
}
compartments<-t(sapply(1:1000, k))

Z1.2<-vectorfun(compartments[,1])
Z2.2<-vectorfun(compartments[,2])

#Let's check the accuracy:
z<-Z1.2==Z1
length(which(z==FALSE))/2

z<-Z2.2==Z2
length(which(z==FALSE))/2
#99.4 and 97.8%, respectively

#Now for computing the Z matrix as probabilities rather than in binary form...

#optimization section: using optim?
#Z1.2*psi1+Z2.2*psi2=phipick
l<-function(psi1){
  phimatrix<-psi1*Z1.2+(1-psi1)*Z2.2
  sum((phimatrix-phipick)^2)
}

optim(c(.75), fn=l, method = "BFGS")
