library(quantmod)
setwd('F:\\H&N')
source("00fracture.R")
#cnv <- na.omit(read.table(file='C:\\Users\\Mark\\Downloads\\hnscc.copy_number.tsv'))
#load('cnv2.rda')
load('hn2.rda')
#cnv <- imputedCNV2
locns <- read.table(file='refGene.hg19.genes.bed', sep="\t", header=FALSE) 
colnames(locns) <- c("Chr", "Start", "End", "Symbol", "Blank", "Strand") 
rownames(locns) <- as.character(locns$Symbol) 
locns$Center <- (locns$Start + locns$End)/2 
locns$Chr <- factor(locns$Chr, levels=paste("chr", c(1:22, "X", "Y"), sep=''))
oo <- order(locns$Chr, locns$Center)
imputedCNV2 <- imputedCNV2[oo,]
gl <- GenomicLocation(locns)
pg <- setupGrid()

#cnv <- na.omit(cnv)
cidlist <- list('GS1445', 'GS1461')
nulist <- list(.5,.5)
#psilist <- c(.73, .53)
#psis <- unlist(psilist)
#threshold <- .05

mixfun <- function(cidlist, nulist, threshold=.05, cutoff=.05){
  nuvec <- unlist(nulist)
  cidvec <- unlist(cidlist)
# psivec <- unlist(psilist)
  jset <- sapply(1:length(cidvec),function(i){which(colnames(imputedCNV2)==cidvec[i])})
  foo <- function(j){
    frac <- Fracture(imputedCNV2[,j], gl, pg, 
         theta=0.9, sigma=0.04, LABEL=j)
    frac@fit$par
  }
  pars <- t(sapply(jset, foo))
  psis <- pars[,2]
  truepsis <- psis*nuvec
  comps <- c(-2, -1, 0, 1, 2, 3, 4)
  datalist <- lapply(1:length(cidvec), function(i){imputedCNV2[,which(colnames(imputedCNV2)==cidvec[i])]})
  midfun <- function(data){
   #den <- density(data)
   #vals <- den$x[findValleys(den$y)]
   #midvals <- c(max(vals[vals < 2]), min(vals[vals > 2]))
   #median(data[data > midvals[1] & data < midvals[2]])
   median(data)
  }
  #mids <- sapply(1:length(datalist), function(i){midfun(datalist[[i]])})
  mids <- pars[,1]
  newlist <- NULL
  for(i in 1:length(datalist)){
    newlist[i] <- list(unlist(datalist[[i]]) + max(mids) - mids[i])
  }
  mixture <- do.call('+', lapply(1:length(cidvec), function(i){nuvec[i]*newlist[[i]]}))
  mid <- midfun(mixture)
  zfun <- function(i){
    lapply(1:length(comps), 
      function(j){which(newlist[[i]] < (mids[i]+psis[i]*comps[j] + threshold) & 
            newlist[[i]] > (mids[i]+psis[i]*comps[j] - threshold))})
  }
  complist <- lapply(1:length(cidlist), zfun)
#finding intersections...
mat <- matrix(rep(NA, 4*length(comps)^2), nrow=length(comps)^2, ncol=4)
  for(i in 1:length(comps)){
    for(j in 1:length(comps)){
      x <- truepsis[1]*comps[i] + truepsis[2]*comps[j] + mid
      y <- length(intersect(complist[[1]][[i]], complist[[2]][[j]]))
      row <- (i-1)*length(comps)+j
      mat[row,1] <- comps[i]
      mat[row,2] <- comps[j]
      mat[row,3] <- x
      mat[row,4] <- y
    }
  }
colnames(mat) <- c('cnv1', 'cnv2', 'peak location', 'number of genes')
sigma <- sd(mixture)
mu <- abs(pnorm(mat[i,3]+threshold) - pnorm(mat[i,3]-threshold))
probs <- unname(sapply(1:nrow(mat),function(i){dbinom(mat[i,4], nrow(imputedCNV2),
                abs(pnorm(mat[i,3]+threshold) - pnorm(mat[i,3]-threshold)))}))
keep <- intersect(which(probs<cutoff), which(mat[,4]>=3))
gridlines <- as.vector(mat[keep,3])
plot(density(mixture), xlim=c(0,4), main=paste(cidvec[1],": ",nuvec[1],"; ",cidvec[2],": ",nuvec[2],sep=''))
abline(v=gridlines, col='red')
list(mat, mixture, truepsis)
}

#temp <- mixfun(list('GS1445','GS1461'),list(.5,.5),)


########Appendix#######
#a <- which(colnames(cnv)=='GS1445')
#data <- cnv[,a]
#den <- density(data)
#vals <- den$x[findValleys(den$y)]
#midvals <- c(max(vals[vals < 2]), min(vals[vals > 2]))
#mid <- median(data[data > midvals[1] & data < midvals[2]])
#plot(density(cnv[,a]),xlim=c(1,3),main=colnames(cnv)[a])
#abline(v=mid+.73*c(-1,1, -2, 2),col='red')

#b <- which(colnames(cnv)=='GS1461')
#data <- cnv[,b]
#den <- density(data)
#vals <- den$x[findValleys(den$y)]
#midvals <- c(max(vals[vals < 2]), min(vals[vals > 2]))
#mid <- median(data[data > midvals[1] & data < midvals[2]])
#plot(density(cnv[,b]),xlim=c(1,3),main=colnames(cnv)[b])
#abline(v=mid+.53*c(-1,1, -2, 2),col='red')
#psi.b <- .53
