rm(list=ls())
source("F:\\Methods2\\CloneFinderFolder\\00-generics.r")
source("F:\\Methods2\\CloneFinderFolder\\01-cloneMaker.r")
source("F:\\Methods2\\CloneFinderFolder\\02-prefit.r")
source("F:\\Methods2\\CloneFinderFolder\\03-psi.r")
source("F:\\Methods2\\CloneFinderFolder\\SNPSim5.r")
#source("F:\\Methods2\\CloneFinderFolder\\working.r")
library(gtools)
library(CloneFinder)
library(stats)
library(igraph)

sigma2.lrr <- .05
sigma2.baf <- .001
xy <- data.frame(x = log10(c(2, 2, 1, 3, 4)/2),
                 y = c(1/2, 0, 0, 1/3, 1/4))

ones <- list(list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1),
             list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1)
             , list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1), list(1))
twos <- list(list(.1, .9), list(.2, .8), list(.25, .75), list(.3,.7), list(.4, .6), list(.45, .55), list(.5, .5),
             list(.08, .92), list(.61, .39), list(.67, .33), list(.53, .47), list(.38, .62), list(.59, .41),
             list(.06, .94), list(.05, .95), list(.035, .965), list(.17, .83), list(.48, .52)
             , list(.18, .82), list(.14, .86), list(.11, .89), list(.23, .77), list(.26, .74), list(.21, .79)
             , list(.44, .56), list(.41, .59), list(.33, .67), list(.31, .69), list(.37, .63), list(.49, .51), list(.45, .55)
             , list(.04, .96), list(.07, .93), list(.9, .11), list(.2, .8), list(.2, .8))
threes <- list(list(.1,.2,.7),list(.1,.3,.6),list(.2,.3,.5),list(.25,.3,.45),list(.28,.34,.38), list(.4, .33, .27),
               list(.1, .14, .77), list(.44, .39, .17), list(.58, .23, .19), list(.07, .33, .6), list(.05, .08, .87)
               ,list(.12,.29,.59),list(.14,.28,.58),list(.19,.25,.64),list(.04,.44,.52),list(.07,.1,.83),list(.16,.35,.49)
               ,list(.32,.27,.41),list(.64,.2,.16),list(.19,.23,.68),list(.05,.41,.54),list(.07,.17,.76),list(.16,.39,.45))
fours <- list(list(.1, .15, .23, .52), list(.07, .15, .37, .41), list(.1, .27, .3, .33), list(.07, .13, .165, .635)
              , list(.12, .18, .26, .44), list(.18, .15, .23, .44), list(.22, .13, .26, .39), list(.06, .21, .29, .44)
              , list(.10, .19, .37, .34), list(.32, .08, .21, .39))
psiList <- c(ones, twos, threes, fours)

###
setwd('D:\\Mark\\Simulations7')
stuff3 <- setGen(psiList)
save(stuff3, file="stuff3.rda")
load('stuff3.rda')
sampleSet <- stuff3[[1]]
conmat <- stuff3[[2]]
confun <- function(data){
  nsegs <- nrow(conmat)
  mat <- matrix(rep(NA, nsegs^2), nrow=nsegs, ncol=nsegs)
  cloneset <- unlist(lapply(1:length(sampleSet), function(i){temp <- data[[i]][[1]][[2]]}), recursive=FALSE)
  nonzeros <- which(sapply(1:nsegs, function(k){length(which(sapply(1:length(cloneset), 
        function(x){cloneset[[x]][k,1]!=1})))>0}))
  for(i in 1:nsegs){
    for(j in 1:nsegs){
      num <- length(which(sapply(1:length(cloneset), function(x){cloneset[[x]][i,1]!=1 && cloneset[[x]][j,1]!=1})))
      denom <- length(which(sapply(1:length(cloneset), function(x){cloneset[[x]][i,1]!=1})))
      mat[i,j] <- num/denom
    }
  }
  mat <- mat[c(nonzeros),c(nonzeros)]
  colnames(mat) <- nonzeros
  rownames(mat) <- nonzeros
  mat
}
conmat.real <- confun(sampleSet) 
#which(sampleSet[[5]][[1]][[2]][[1]][,1]==0)
heatmap(conmat, main="Co-occurrence parameters")
heatmap(conmat.real, main="Actual co-occurrence")

analyze <- function(n){
  obj <- sampleSet[[n]] 
  #obj <- genomegen(30, threes[[5]], list(.3, .3, .3, .1), .5, 600000, 10000, .05, .001, seed=1)
  markers <- obj[[2]][[2]]
  dataset <- obj[[2]][[1]]
  wholedata <- obj[[3]][[1]]
  colnames(dataset) <- c("x", "y")
  dataset <- as.data.frame(dataset)
  #dataset <- data.frame("x"=obj[[2]][[1]], "y"=obj[[2]][[2]])
  truepsis <- rev(unlist(obj[[1]][[1]]))
  marker.data <- Reduce("+", lapply(1:length(truepsis), function(i){truepsis[i]*obj[[3]][[i]]}))
  ZList <- obj[[1]][[2]]
  #bafpts <- obj[[3]]
  #lrrpts <- obj[[4]]
  
  compModel <- CompartmentModel(markers, xy, sigma2.lrr)
  
  # good guess at psi-vector
  pcm <- PrefitCloneModel(dataset, compModel)
  upd <- updatePhiVectors(pcm, compModel)
  mfun <- function(m){
    set.seed(n+3)
    # temp <- upd@phipick
    # for(i in 1:nrow(temp)){
    #   inds <- which(temp[i,] %in% sort(temp[i,], decreasing=TRUE)[1:n])
    #   exc <- which(!1:5 %in% inds)
    #   vals <- temp[i,inds]
    #   temp[i,inds] <- vals/sum(vals)
    #   temp[i,exc] <- 0
    # }
    #upd@phipick <- temp
    estpsi <- guessPsi(upd, m) # 3 is number of clones we are trying to fit
    #final <- runEMalg(estpsi, dataset, compModel, upd@phipick)
    #not sure why upd@phipick is there?
    final <- runEMalg(estpsi, dataset, compModel)
    list(final, upd@phipick)
  }
  results <- sapply(1:6, function(z){mfun(z)[[1]]})
  truth <- list(truepsis, ZList, markers)
  output <- list(results, truth, dataset, lapply(1:6, function(z){mfun(z)[[2]]}), wholedata)
#  assign(paste("sim.", n, sep=""),output)
  #save(list=paste("sim.", n, sep=""), file=paste("D:\\Mark\\Simulations5\\sim.", n, ".rda", sep=""))
#  save(output, file=paste("D:\\Mark\\Simulations5\\sim.", n, ".rda", sep=""))
 output
}
res3 <- lapply(1:length(psiList), analyze)
load('results-11-6-15.rda')
save(res3, file="results-11-6-15.rda")
res <- res3
#start: 1:43

###Co-occurrence analysis:
nrows <- nrow(res[[1]][[1]][2,][[1]])
nsamples <- length(res)

picked <- unlist(lapply(1:length(res), function(i){min(which(sapply(1:6, function(j){min(unlist(
  res[[i]][[1]][1,j]))})<.04))}))-1
Z.temp <- lapply(1:length(res), function(i){res[[i]][[1]][2,picked[i]]})
priorfun <- function(q){
  inner <- function(n){
   Zs <- res[[n]][[1]][2,picked[n]]$Zmats
   clonefun <- function(i){
     Zs[q,1,i]!=1
   }
   outp <- sapply(1:picked[n], clonefun)
  }
 bin <- unlist(sapply(1:nsamples, inner))
 length(which(bin))/length(bin)
}
priors <- sapply(1:nrows, priorfun)
sample <- as.vector(res[[1]][[4]][[1]])
sd0 <- sd(sample[sample>.5])
mat <- matrix(rep(NA, nrows^2), nrow=nrows, ncol=nrows)
infer <- function(dataset, s2=sd0){
  inner <- function(x){
    innerer <- function(y){
      innerest <- function(n){
       Zs <- res[[n]][[1]][2,picked[n]]$Zmats
        clonefun <- function(i){
         x.norm <- Zs[x,1,i] 
         y.norm <- Zs[y,1,i]
         status <- NULL
         if(x.norm==0 & y.norm==1){
           status <- "x"
         }else if(x.norm==1 & y.norm==0){
           status <- "y"
         }else if(x.norm==1 & y.norm==1){
           status <- "neither"
         }else if(x.norm==0 & y.norm==0){
           status <- "both"
         }
         status
        }
       sapply(1:picked[n], clonefun)
      }
      tab <- table(unlist(sapply(1:nsamples, innerest)))
      df <- as.data.frame(t(data.frame(as.vector(tab))))
      colnames(df) <- names(tab)
      rownames(df) <- 'values'
      prior <- priors[x]*priors[y]
      if(length(df$x)==0){
        df$x <- 0
      }
      if(length(df$both)==0){
        df$both <- 0
      }
      val <- as.numeric(df$both)/as.numeric(df$both+df$x)
      approx <- function(par){
        dnorm(val, par, s2)*dnorm(val, prior, s2)
      }
      if(prior==0 | length(val)==0){
        mat[x,y] <- NA
      }else{
        parvec <- seq(from=0, to=1, length=101)
        probs <- sapply(parvec, approx)/sum(sapply(parvec, approx))
        EV <- sum(probs*parvec)
        #mat[x,y] <- EV
        mat[x,y] <- val
      }
      val
    }
    sapply(1:nrows, innerer)
  }
  sapply(1:nrows, inner)
}
conmat.inf <- infer(res)
nonzeros <- which(priors>0)
submat.inf <- conmat.inf[c(nonzeros),]
submat.inf <- submat.inf[,c(nonzeros)]
submat.true <- conmat[c(nonzeros),]
submat.true <- submat.true[,c(nonzeros)]

colnames(submat.inf) <- nonzeros
colnames(submat.true) <- nonzeros
rownames(submat.inf) <- nonzeros
rownames(submat.true) <- nonzeros

#compare posts with conmat using heatmaps; also with the actual co-ocurrence freqs
dm.true <- data.matrix(submat.true)
dm.inf <- data.matrix(submat.inf)
dm.real <- data.matrix(conmat.real)
dm.true.sub <- data.matrix(conmat[c(as.numeric(rownames(conmat.real))),c(as.numeric(rownames(conmat.real)))])
heatmap(dm.true)
heatmap(dm.inf)
plot(as.vector(submat.true), as.vector(submat.inf))

heatmap(dm.true.sub)
heatmap(dm.real)

##################Tree Likelihoods########
n <- 81
phidata <- res[[n]][[4]][5][[1]]
binmat <- phidata > .04 & phidata < .96
nonzero <- which(binmat) 
cols <- ceiling(nonzero/nrow(phidata))
rows <- nonzero - (cols-1)*nrow(phidata)
mat <- conmat.inf*nrow(phidata)
cooc <- matrix(rep(NA, length(cols)^2), nrow=length(cols), ncol=length(cols))
for(i in 1:length(cols)){
  for(j in 1:length(cols)){
    if(cols[i]>1 & cols[j]>1){
      cooc[i,j] <- conmat.inf[i,j]
    }else if(cols[i]==1 & cols[j]>1){
      cooc[i,j] <-  priors[j] - conmat.inf[i,j]*priors[j]
    }else if(cols[i]>1 & cols[j]==1){
      cooc[i,j] <- priors[i] - conmat.inf[i,j]*priors[j]
    }else{
      cooc[i,j] <- priors[i] + priors[j] - conmat.inf[i,j]*priors[j] 
    }
  }
}
phi <- sort(phidata[nonzero], decreasing=TRUE)

N <- nrow(cooc)
L <- factorial(N)
ilists <- list(NULL)
for(i in 1:L){
  ilists[i] <- list(rep(0, N))
}

pool <- 0:(N-1)
i <- 2
stems <- list(rep(0, N))
ilist <- stems
while(i <= N) {
  open <- pool[pool<i]
  leaves <- list(NULL)
  for(j in 1:length(stems)){
    for(k in open){
      leaf <- stems[[j]]
      leaf[i] <- k
      leaves <- c(leaves, list(leaf))
    }
  }
  leaves <- leaves[2:length(leaves)]
  stems <- leaves
  i <- i+1
}
ilists <- stems
N <- N-1

prods <- NULL
probs <- NULL
probs2 <- NULL
sigma <- .1
sigma2 <- .1
for(i in 1:L){
  ivec <- ilists[[i]]
  dmat <- matrix(rep(NA, 2*(N+1)), nrow=2, ncol=N+1)
  dmat[,1] <- c(1, phi[1])
  for(k in (N+1):2){
    dmat[,k] <- c(phi[ivec[k-1]], phi[k])
  }
  psi <- c(dmat[1,] - dmat[2,], phi[length(phi)])
  prods[i] <- prod(cooc[c(1:N)+1,c(ivec)])
  ivec <- c(0, ivec)
  phivec <- c(1, phi[1:(length(phi)-1)])
  res <- NULL
  for(j in N:0){
    res[j] <- phivec[j+1] - sum(phi[which(ivec==j)])
  }
  less <- which(res<0)
  res2 <- res[less]/phi[less]
  if(length(res2)>0){
    probs[i] <- prod(pnorm(res2, 0, sigma))
  }else{
    probs[i] <- .5
  }
  probs2[i] <- dnorm(1-sum(psi), 0, sigma2)
}
imax <- ilists[[which.max(prods*probs*probs2)]]

dmat <- matrix(rep(NA, 2*(N+1)), nrow=2, ncol=N+1)
dmat[,1] <- c(1, phi[1])
for(i in (N+1):2){
  dmat[,i] <- c(phi[imax[i-1]], phi[i])
}

psi <- c(dmat[1,] - dmat[2,], phi[length(phi)])
tab <- rbind("parent"=c(NA, imax), psi)
tab <- rbind(tab, "rows"=c(1, rows), "cols"=c(NA, cols))

mat <- matrix(rep(0, nrows*5), ncol=5, nrow=nrows)
mat[,1] <- rep(1, nrows)
startmat <- mat
mats <- list(startmat)
i <- 1
while(i<=N+2) {
  if(i>1){ 
    Z <- mats[[tab[1,i]+1]]
    Z[as.integer(tab[3,i]),] <- rep(0, 5)
    Z[as.integer(tab[3,i]),as.integer(tab[4,i])] <- 1
  }
 if(i>1){
   mats <- c(mats, list(Z))
 }
 i <- i+1
}
which()
############################################################################################

check <- function(k){
  output <- set[[k]]
  sim <- output
  truepsis <- sort(sim[[2]][[1]], decreasing=TRUE)
  psiset <- output[[1]][1,]
  nclone <- min(which(sapply(1:length(psiset), function(i){min(psiset[[i]])<.05})))-1
  phimat <- output[[3]]
  trueZ <- rev(sim[[2]][[2]])[[1]]
  estZ <- sim[[1]][2+3*(length(truepsis)-1)][[1]]
  estZ <- estZ[,,1]
  data <- output[[4]]
  falses <- which(sapply(1:nrow(estZ), function(i){length(which(estZ[i,]==trueZ[i,]))!=5}))
  plot(sapply(1:nrow(phimat), function(i){max(phimat[i,])}))
  #c("true"=length(truepsis), "guess"=nclone)
}

res <- t(sapply(1:length(psiList), check))
length(which(res[,1]==res[,2]))/nrow(res)


###Compare results with the truth:
foo <- function(k){
  output <- set[[k]]
  #  truepsis <- sort(sim[[2]][[1]], decreasing=TRUE)
  #  nclone <- length(truepsis)
  #  estpsis <- sim[[1]][[1+3*(nclone-1)]]
  #  psidiff <- sum((sort(truepsis) - sort(estpsis))^2)
  #  estZ <- sim[[1]][2+3*(nclone-1)][[1]]
  #  trueZ <- rev(sim[[2]][[2]])
  trueZ <- output[[2]][[2]]
  nclone <- length(trueZ)
  estZ <- output[[1]][2,][[nclone]]
  compare <- function(i){
    bin <- trueZ[[i]]==estZ[,,i]
    bin2 <- sapply(1:nrow(estZ), function(j){length(which(bin[j,]==FALSE))>0})
    1-(length(which(bin2==TRUE))/nrow(estZ))
  }
  com <- sapply(1:nclone, compare)
  falses <- lapply(1:nclone, function(n){which(sapply(1:nrow(estZ), function(i){
    length(which(trueZ[[n]][i,]==estZ[i,,n]))!=5}))})
  #  list(t(data.frame("true psis"=truepsis, "inferred"=estpsis)), com)
  com
}
coms <- lapply(1:length(set), foo)

###
#Assuming objects 'res3' and 'stuff3
real <- stuff3[[1]]
thetas <- seq(from=.01, to=.9, length=90)
thetatest <- function(theta){
 trueNs <- sapply(1:length(real),function(i){length(real[[i]][[1]][[1]])})
 inferredNs <- sapply(1:length(res3),function(i){
   which.max(unlist(res3[[i]][[1]][3,])*dexp(1:6,theta))})
 list(length(which(trueNs==inferredNs)), data.frame('true'=trueNs,'inferred'=inferredNs))
}
L <- lapply(thetas, thetatest)
i.best <- which.max(sapply(1:length(L),function(i){L[[i]][[1]]}))
tab <- L[[i.best]][[2]]
best <- thetas[i.best]
acc <- length(which(tab[,1]==tab[,2]))/nrow(tab)
accs <- sapply(1:length(thetas),function(i){
  length(which(L[[i]][[2]][,1]==L[[i]][[2]][,2]))/nrow(L[[i]][[2]])})
plot(thetas, accs,type='l')
