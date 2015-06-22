set.seed(1)
bases<- c("A", "T", "C", "G")
l<- 50000
transitionRate <- 1/20
transversionRate <- 1/40
#l is sequence length
seqGen<- function(n){
  sample(bases, n, replace=FALSE)
}
#note: base content can be modified with prob option (e.g. prob=c(.1, .1, .4, .4))
seq<- t(sapply(rep(4, l), seqGen))

#intensity matrix
v <- c(rbeta(l, 20, 1), rbeta(3*l, 1, 20))
m <- matrix(v, nrow=l, ncol=4)
sumfun<- function(n){
  mat<- matrix(nrow=l, ncol=4)
  mat[n,] <- m[n,]/sum(m[n,])
}
m2<- t(sapply(1:l, sumfun))

#SNPs and mutations:
estBetaParams <- function(mu, s2) {
  temp <- mu*(1-mu)/s2 - 1
  alpha <- mu*temp
  beta <- (1 - mu)*temp
  c(alpha = alpha, beta = beta)
}

CGs <- which(seq[,1] %in% "C" & seq[,2] %in% "G")
GCs <- which(seq[,1] %in% "G" & seq[,2] %in% "C")
TAs <- which(seq[,1] %in% "T" & seq[,2] %in% "A")
ATs <- which(seq[,1] %in% "A" & seq[,2] %in% "T")
transit <- sort(c(CGs, GCs, TAs, ATs))
transvers <- setdiff(1:l, transit)
ratio <- length(transit)/length(transvers)
sample1<- sample(transit, round(transitionRate*length(transit)), replace=FALSE)
sample2<- sample(transvers, round(transversionRate*length(transvers)), replace=FALSE)
SNPindex <- c(sample1, sample2)

SNPGen<- function(n, s=.0001){
  mn<- (sum(m2[n,1])+sum(m2[n,2:4])/3)/2
  v <- estBetaParams(mn, s)
  vals <- rbeta(2, v[1], v[2])
  m2[n,1]<- vals[1]
  m2[n,2]<- vals[2]
  m2[n,]<- m2[n,]/sum(m2[n,])
}

SNPvals <- t(sapply(SNPindex, SNPGen))
m2[SNPindex,] <- SNPvals

#Signal adjustment; using lognormal distribution
w<- rlnorm(l, 5, .1)
signal <- m2*w

#Let's consider a more complex error structure:
allele1<- seq[,1]
allele2<- allele1
allele2[SNPindex] <- seq[,2][SNPindex]
newseqs <- cbind(allele1, allele2)
seqdf <- data.frame("Allele1"=allele1, "Allele2"=allele2)

#CNVs, clones, and both alleles having alt variants:
histclones <- 5
nclones <- 3
elim <- sample(1:histclones, histclones-nclones)
keep <- which(!1:histclones %in% elim)
nclones <- 1
psis <- c(1)
homfreq <- .4
hetfreq <- .4
Gainfreq <- .1
Lossfreq <- .1
#hom=1, het=0
CNVmed <- l/4
CNVmin <- l/25
clonegen <- function(baseclone){
  start <- baseclone
  ind <- round(runif(1, 1, l))
  allel <- sample(c(0, 1, 2, 3), 1, prob=c(homfreq, hetfreq, Gainfreq, Lossfreq))
  basecols <- ncol(baseclone)
  if(allel==0){
    newbase <- sample(bases[bases != unname(as.matrix(seqdf)[ind,])], 1)
    start[ind,] <- newbase
  }else if(allel==1){
    toRep <- sample(c(1, 2), 1)
    newbase <- sample(bases[bases != unname(as.matrix(seqdf)[ind,toRep])], 1)
    start[ind,toRep] <- newbase
  }else if(allel==3){
    len <- abs(rnorm(1, CNVmed, 10))+CNVmin
    toChange <- intersect((ind-round(len/1)):(ind+round(len/1)), 1:l)
    pick <- sample(1:2, 1)
    start[,basecols+1] <- as.character(rep("NA", l))
    start[toChange,basecols+1] <- as.character(start[toChange,pick])
  }else{
    len <- abs(rnorm(1, CNVmed, 10))+CNVmin
    toChange <- intersect((ind-round(len/1)):(ind+round(len/1)), 1:l)
    pick <- sample(1:2, 1)
    start[which(!1:l %in% toChange),pick] <- 0
  }
  start
}
set.seed(12)
n <- 1
clones <- list(seqdf)
while(n < histclones) {
  pick <- sample(1:n, 1)
  newclone <- clonegen(clones[[pick]])
  clones <- unlist(list(clones, list(newclone)), recursive=FALSE)
  n <- n+1
}
clones <- clones[keep]

variantfind <- function(i){
  innerfun <- function(n){
    row <- clones[[n]][i,][clones[[n]][i,]!="NA"]
    length(row)==2 & length(unique(row))==1
  }
  length(which(sapply(1:length(clones), innerfun)==TRUE))==length(clones)
}
allvariants <- which(sapply(1:nrow(seqdf), variantfind)!=TRUE)
vtypes <- c("Homozygous SNP", "Heterozygous SNP", "Subclonal Mutation", "CNV")
variant.type <- 
  
outerfun <- function(clone){
seqdf <- clones[[clone]]
signalFun <- function(i){
  b<- 100
  #b is the base signal parameter; must be at least as large as read length.
  m <- 100
  #m is a scaling parameter that relates to the mean signal intensity.
  d <- seqdf[1:i,]
  n <- 80
  r<- .2
  #r is the base noise level
  s<-.08
  #n is a scaling parameter
  signalA <- m*sum(((b - (i-which(d[,1] %in% "A")))/b)^n)
  signalT <- m*sum(((b - (i-which(d[,1] %in% "T")))/b)^n)
  signalC <- m*sum(((b - (i-which(d[,1] %in% "C")))/b)^n)
  signalG <- m*sum(((b - (i-which(d[,1] %in% "G")))/b)^n)
  #"base" noise is defined by the term: rlnorm(l, m, s)
  ObsA<- rlnorm(1, log(signalA), s) + rlnorm(1, log(r), s)
  ObsT<- rlnorm(1, log(signalT), s) + rlnorm(1, log(r), s)
  ObsC<- rlnorm(1, log(signalC), s) + rlnorm(1, log(r), s)
  ObsG<- rlnorm(1, log(signalG), s) + rlnorm(1, log(r), s)
  c(ObsA, ObsT, ObsC, ObsG)
}
ObservedMat <- t(rbind(sapply(1:80, signalFun)))
dfObserved <- data.frame("A"=ObservedMat[,1], "T"=ObservedMat[,2], "C"=ObservedMat[,3], "G"=ObservedMat[,4])
psis[clone]*dfObserved
}
if(nclones <- 1){
  psis <- 1
}
netsignal <- round(Reduce("+", lapply(1:nclones, outerfun)))

ObsFun<- function(n){
  max(netsignal[n,])
}
maxes <- sapply(1:80, ObsFun)

plot(1:80, maxes)
lines(netsignal[,1], col="red")
lines(netsignal[,2], col="blue")
lines(netsignal[,3], col="green")
lines(netsignal[,4], col="yellow")

