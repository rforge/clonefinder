 #source("F:\\Methods2\\betamix-emalg.R")
#l <- 600000
#h <- 20
#mu.l <- 10000
#psilist <- ones[[1]]
#psilist <- truepsis
#lrr.s2 <- .05
#baf.s2 <- .001
#probs <- list(.3, .3, .3, .1)
#snprate <- .5
#seed <- 1

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

setGen <- function(psiList, lrr.sigma=.05^.5, baf.sigma=.001^.5, nsegs=52, l=600000, mu.l=10000, q=25, probs=list(.3, .3, .3, .1)){
abset <- 1:(nsegs-21)
cbrk <- runif(22, 1, 2)
cbrk <- sort(round((cbrk/sum(cbrk))*l))
cbrk <- sapply(1:22, function(i){sum(cbrk[1:i])})
cbrk[22] <- l
starts <- c(1, sort(round(runif(nsegs-22, 0, l))))
lens <- round(rlnorm(nsegs-21, log(mu.l), log(2.5)))

x <- log10(c(2, 1, 3, 4)/2)
y <- c(0, 0, 1/3, 1/4)

#co-occurrence values, as follows: 
prior.mu <- .1
prior.var <- 50
alpha <- abs(unlist(estBetaParams(prior.mu, prior.var)[1]))
beta <- abs(unlist(estBetaParams(prior.mu, prior.var)[2]))
conmat <- matrix(rbeta(nsegs^2, alpha, beta), nrow=nsegs, ncol=nsegs)
conmat[conmat>1] <- 1/conmat[conmat>1]
conmat <- abs(conmat)
for(i in 1:nsegs){
 conmat[i,i] <- 1
}

genomegen <- function(n, snprate=.5, lrr.s2=lrr.sigma, baf.s2=baf.sigma, seed=1){
psivec <- unlist(psiList[[n]])
h <- sample(length(psivec):q, 1, prob=1/length(psivec):q)
CNVs <- sample(1:4, h, prob=probs, replace=TRUE)
probs <- unlist(probs)
psis <- sort(unlist(psivec))
set.seed(seed)

sets <- list()
i <- 0
while(i < h) {
  set.seed(n*h+i)
  pool <- which(!1:nsegs %in% unique(unlist(sets)))
  ind <- sample(1:length(sets), 1)
  base <- unlist(sets[ind])
  probmat <- conmat[pool,base]
  if(length(base)>1){
    probvec <- apply(probmat, 1, prod) 
  }else{
    probvec <- probmat
  }
  probvec <- probvec/sum(probvec)
  if(length(base)>0){
    pick <- sample(pool, 1, prob=probvec)
  }else{
    pick <- sample(pool, 1)
  }
  newclone <- list(base, pick)
  sets <- unlist(list(sets, list(newclone)), recursive=FALSE)
  i <- length(sets)
}

picked <- sample(1:h, length(psis))
ends <- sapply(1:length(starts), function(j){min(starts[j]+lens[j], min(cbrk[cbrk-starts[j]>0]), starts[j+1]-1)})
ends[nsegs-21] <- l
#seg.start <- sort(c(starts[unique(unlist(sets[picked]))], c(1,cbrk[1:21]+1)))
seg.start <- sort(c(starts, cbrk[1:21]+1))
#seg.end <- sort(c(ends[unique(unlist(sets[picked]))], c(cbrk)))
seg.end <- sort(c(starts[2:length(starts)]-1, c(cbrk)))
segs <- t(rbind(seg.start, seg.end))

foo <- function(i, alpha=.1){
  xdata <- rep("NA", l)
  ydata <- rep("NA", l)
  for(j in unique(unlist(sets[[picked[i]]]))){
  #end1 <- starts[j]+lens[j]
  #end2 <- min(cbrk[cbrk-starts[j]>0])
  #end3 <- starts[j+1]-1
  #endpt <- min(end1, end2, end3)
  endpt <- seg.end[j]
  absegs <- sort(unique(unlist(sets)))
  q <- which(absegs %in% j)
  xdata[seg.start[j]:endpt] <- rnorm(length(seg.start[j]:endpt), x[CNVs[q]], lrr.s2)
  ydata[seg.start[j]:endpt] <- rnorm(length(seg.start[j]:endpt), y[CNVs[q]], baf.s2)
  }
  xdata[is.na(xdata)] <- rnorm(length(xdata[is.na(xdata)]), 0, lrr.s2)
  xdata[xdata=="NA"] <- rnorm(length(xdata[xdata=="NA"]), 0, lrr.s2)
  ydata[is.na(ydata)] <- rnorm(length(ydata[is.na(ydata)]), .5, baf.s2)
  ydata[ydata=="NA"] <- rnorm(length(ydata[ydata=="NA"]), .5, baf.s2)
  ydata <- as.numeric(ydata)
  xdata <- as.numeric(xdata)
  ydata[ydata<0] <- -ydata[ydata<0]
  ydata[ydata>1] <- 1/ydata[ydata>1]
  noise <- sample(1:l, round(l*(1-snprate)))
  ydata[noise] <- rbeta(round(l*(1-snprate)), alpha, 100)
  minus <- runif(round(l/2), 0, l)
  ydata[minus] <- 1-ydata[minus]
  dat<- t(rbind(xdata, ydata))
  dat
}
data <- lapply(1:length(psis), foo)

foo2 <- function(n){
  data <- data[[n]]
  innerfun1 <- function(i){
    median(data[seg.start[i]:seg.end[i],1])
  }
  lrr <- sapply(1:length(seg.start), innerfun1)
  baf <- data[,2]
  baf[baf>.5] <- 1- baf[baf>.5]
  bar <- .1
  innerfun2 <- function(i){
    stuff <- baf[seg.start[i]:seg.end[i]]
#    betamix(stuff[stuff>.1])
    if((length(stuff[stuff>.1])/length(stuff[stuff<.1]))>.25){
      med <- median(stuff[stuff>.1])
    }else{
      med <- median(stuff)
    }
    med
  }
  baf <- sapply(1:length(seg.start), innerfun2)
  t(rbind(lrr, baf))
}
segdata <- lapply(1:length(psis),foo2)
segmentdata <- Reduce("+", lapply(1:length(psis), function(i){segdata[[i]]*psis[i]}))

Zfun <- function(j){
Z <- matrix(rep(0, 5*nrow(segs)), nrow=nrow(segs), ncol=5)
Z[,1] <- 1
inds <- which(seg.start %in% starts[unlist(sets[picked[j]])])
Z[inds,1] <- 0
cols <- CNVs[sort(unlist(sets[picked[j]]))]+1
for(i in 1:length(inds)){
  Z[inds[i],cols[i]] <- 1
}
Z
}
Zs <- lapply(1:length(psis), Zfun)

markers <- seg.end-seg.start+1
list(list(psis, Zs), list(segmentdata, markers), data)
}
sampleSet <- lapply(1:length(psiList), genomegen, snprate=.5, lrr.s2=.05, 
                    baf.s2=.001, seed=1)
list(sampleSet, conmat)
}

infer <- function(table){
  
}
