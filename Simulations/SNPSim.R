set.seed(1)
hclones <- 6
nclones <- 3
psis <- c(.5, .3, .2)
#l is the number of markers
l<- 30000
#The two 'alleles' are sets of markers along the paternal or maternal chromosome sets
allele.1 <- rep("A", l)
allele.2 <- rep("A", l)
#snprate is the rate of snps across the marker set
snprate <- .5
bset1 <- round(runif(round(l*snprate), 1, l))
bset2 <- round(runif(round(l*snprate), 1, l))
allele.1[bset1] <- "B"
allele.2[bset2] <- "B"
basesnp <- data.frame(allele.1, allele.2)
mu <- 5000

lossfun <- function(seq){
  lvec <- 1:nrow(seq)
  len <- round(rlnorm(1, log(mu), .5))
  ind <- round(runif(1, 1, nrow(seq)-len))
  allele <- sample(1, 1:ncol(seq))
  newseq <- as.matrix(seq)
  newseq[ind:(ind+len),allele] <- rep(0, len+1)
  newseq
}

gainfun <- function(seq){
  lvec <- 1:nrow(seq)
  len <- round(rlnorm(1, log(mu), .5))
  ind <- runif(1, 1, nrow(seq)-len)
  allele <- sample(1, 1:ncol(seq))
  newseq <- as.matrix(seq)
  gain <- newseq[ind:(ind+len),allele]
  newseq <- as.data.frame(newseq)
  newseq[,(ncol(seq)+1)] <- rep(0, nrow(newseq))
  newseq <- as.matrix(newseq)
  newseq[ind:(ind+len),(ncol(seq)+1)] <- gain
  newseq
}

#by default we'll assume losses and gains are equally likely
CNVfun <- function(clone, probvec=c(.5, .5)){
  sam <- sample(1:2, 1, prob=probvec)
  if(sam==1){
    newclone <- lossfun(clone)
  }else{
    newclone <- gainfun(clone)
  }
  newclone
}
clones <- list(basesnp)

while((length(clones)-1) < hclones) {
  i <- sample(1, 1:length(clones))
  baseclone <- clones[[i]]
  newclone <- CNVfun(baseclone)
  clones <- unlist(list(clones, list(newclone)), recursive=FALSE)
}
clones <- clones[2:length(clones)]
extant <- clones[sample(1:length(clones), nclones)]

#Generating theoretical LRR and BAF
BAFfun <- function(clone){
  df <- extant[[clone]]
  innerfun <- function(i){
    row <- unname(df[i,])
    row <- row[row!="0"]
    length(row[row=="B"])/length(row)
  }
  psis[clone]*sapply(1:nrow(df), innerfun)
}

LRRfun <- function(clone){
  df <- extant[[clone]]
  innerfun <- function(i){
    row <- unname(df[i,])
    row <- row[row!="0"]
    log10(length(row)/2)
  }
  psis[clone]*sapply(1:nrow(df), innerfun)
}

BAF <- Reduce("+", lapply(1:nclones, BAFfun))
LRR <- Reduce("+", lapply(1:nclones, LRRfun))

#Introducing noise:
estBetaParams <- function(mu, s2) {
  temp <- mu*(1-mu)/s2 - 1
  alpha <- mu*temp
  beta <- (1 - mu)*temp
  abs(c(alpha = alpha, beta = beta))
}

bafgen <- function(i, s2=.02){
  mu <- BAF[i]
  bafs <- rnorm(1, mu, s2)
  bafs[bafs < 0] <- -bafs[bafs < 0]
  bafs[bafs > 1] <- 1 - (bafs[bafs > 1]-1)
  bafs
}
dat.baf <- sapply(1:length(BAF), bafgen)

lrrgen <- function(i, s2=.02){
  mu <- LRR[i]
  rnorm(1, mu, s2)
}
dat.lrr <- sapply(1:length(LRR), lrrgen)
plot(dat.lrr, cex=.01)
plot(dat.baf, cex=.01)

#Segmentation:
v <- var(dat.lrr)
m <- median(dat.lrr)
breakfun <- function(i, interval=25){
  mdiff <- abs(median(dat.lrr[(i-interval):i]) - median(dat.lrr[i:(i+interval)]))
  sddiff<- abs(sd(dat.lrr[(i-interval):i]) - sd(dat.lrr[i:(i+interval)]))
  mdiff/sddiff
}

ratios <- sapply(51:(l-50), breakfun)
plot(ratios, type="l")
mean(ratios)
threshold <- 3*sd(ratios)
changepoints <- which(ratios>threshold)
startvec <- c(1, changepoints)
endvec <- c(changepoints, l)
yvec <- sapply(1:length(startvec), function(i){median(dat.lrr[startvec[i]:endvec[i]])})
plot(dat.lrr, cex=.01)
segments(x0=startvec, y0=yvec,, x1=endvec, col="green")
markers <- endvec-startvec
med.lrr <- yvec
dat.baf.fold <- dat.baf
dat.baf.fold[dat.baf.fold>.5] <- 1 - dat.baf.fold[dat.baf.fold>.5]
med.baf <- sapply(1:length(startvec), function(i){
  median(dat.baf.fold[startvec[i]:endvec[i]])})
list(med.lrr, med.baf, markers)
