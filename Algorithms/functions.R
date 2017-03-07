library(gtools)
#load('na.array.rda')

logLike <- function(resids, nmarks, sigma0){
  dnorm(resids, 0, sigma0/(nmarks)^.5, log=TRUE)
}
#Length of resids = length of nmarks

psiPrior <- function(psi, alpha){
  log(ddirichlet(psi[c(which(psi>0))], rep(alpha, length(which(psi>0)))))
}

cnPrior <- function(cnmat, theta){
  sum(dgeom(abs(1 - as.vector(cnmat)), theta, log=TRUE))
}

cnPrior.vec <- function(cnvec, theta){
  sum(dgeom(abs(1 - cnvec), theta, log=TRUE))
}

kPrior <- function(k, ktheta){
  dgeom(k, ktheta, log=TRUE)
}

crossover <- function(parents1, parents2, O){
  mats <- lapply(1:length(parents1), function(i){
    par1 <- parents1[[i]]
    par2 <- parents2[[i]]
    breakpoints <- sample(2:ncol(par1), O, replace=FALSE)
    lapply(1:length(breakpoints), function(j){
      cbind(par1[,1:(breakpoints[j]-1)], par2[,breakpoints[j]:ncol(par2)])
    })
  })
  mats <- unlist(mats, recursive=FALSE)
  mats
}
#note: keep in mind alternative method for crossover of sampling rows from 
#each parent rather than just picking a breakpoint.

crossover2 <- function(parents1, parents2, O){
  mats <- lapply(1:length(parents1), function(i){
    par1 <- parents1[[i]]
    par2 <- parents2[[i]]
    n1 <- sample(1:nrow(par1), 1)
    lapply(1:length(1:O), function(j){
      frompar1 <- sample(1:nrow(par1), n1)
      frompar2 <- which(!1:nrow(par1) %in% frompar1)
      output <- matrix(data=NA, nrow=nrow(par1), ncol=ncol(par1))
      for(j in frompar1){
        output[j,] <- par1[j,]
      }
      for(j in frompar2){
        output[j,] <- par2[j,]
      }
      output
    })
  })
  mats <- unlist(mats, recursive=FALSE)
  mats
}

findPsi <- function(datavec, AB){
  AB <- as.matrix(AB)
  M <- t(AB)
  keep <- which(sapply(1:ncol(M), function(i){length(which(M[,i]==0)) < nrow(M)}))
  M <- M[,c(keep)]
  dat <- as.matrix(datavec)
  psi <- as.vector(solve(t(M) %*% M) %*% t(M) %*% dat)
  psi <- abs(psi)
  psi <- psi/sum(psi)
  psi
}
#(M(T)M)^(-1)*M(T)*W; W = c(X, Y) 

convert <- function(log.probs){
  probs.temp <- exp(log.probs)
  if(length(which(probs.temp>0))==length(log.probs)){
    probs <- probs.temp
    probs <- probs/sum(probs)
  }else{
    maxim <- max(log.probs)
    probs <- log.probs - maxim
    probs <- exp(probs)
    minim <- min(probs)
    probs[which(probs==0)] <- minim*(1/sum(probs) + length(which(probs==0))*minim)
    probs <- probs/sum(probs)
    if(length(which(probs==0))==length(probs)){
      probs <- rep(1, length(probs))
      probs <- probs/sum(probs)
    }else if(length(which(probs==1))==1){
      probs[which(probs==0)] <- 1/length(which(probs==0))
      probs <- probs/sum(probs)
    }
  }
  probs
}

#N must be less than or equal to nrow(psipool)
filter <- function(psipool, vecpool, datavec){
  N <- nrow(psipool)
  arr <- vecpool %*% t(psipool)
  matlist <- lapply(1:nrow(psipool), function(i){
    mat <- t(sapply(1:length(datavec), function(j){
      r <- (arr[,i] - datavec[j])^2
      minim <- min(r)
      index <- which.min(r)
      c(minim, index)
    }))
  })
  rsums <- sapply(1:length(matlist), function(i){sum(matlist[[i]][,1])})
  vec.array <- array(data=NA, dim=c(nrow(psipool), length(datavec), ncol(psipool)))
  for(i in 1:nrow(psipool)){
    vec.array[i,,] <- vecpool[c(matlist[[i]][,2]),]
  }
  residmat <- t(sapply(1:length(matlist), function(i){matlist[[i]][,1]}))
  list(psipool, vec.array, residmat)
}

filter2 <- function(psipool, vecpool, datavec){
  for(i in 1:nrow(psipool)){
    
  }
}

resample.grid <- function(base.grid, posts, N){
  resam <- base.grid[c(sample(1:nrow(base.grid), N, replace=TRUE, prob=convert(posts))),]
  mus <- colMeans(resam)
  sums <- sapply(1:nrow(resam), function(i){sum(resam[i,])})
  ns <- sample(sums, N, replace=TRUE)
  ps <- t(sapply(1:length(ns),function(i){mus/ns[i]}))
  t(sapply(1:length(ns),function(i){
    if(ns[i]>0){
      out <- as.vector(rmultinom(1, ns[i], ps[i,]))
    }else{
      out <- rep(0, ncol(resam))
    }
    out
  }))
}

resample.psi <- function(base.psimat, posts, N){
  base.psis <- base.psimat
  for(i in 1:nrow(base.psis)){
    base.psis[i,which(base.psis[i,]<=.02)] <- 0
    base.psis[i,] <- base.psis[i,]/sum(base.psis[i,])
  }
  resam <- base.psis[c(sample(1:nrow(base.psis), N, replace=TRUE, prob=convert(posts))),]
  nonzero <- sapply(1:ncol(resam), function(i){length(which(resam[,i]>0))})
  mus <- sapply(1:ncol(resam),function(i){sum(resam[,i])/nonzero[i]})
  vars <- sapply(1:ncol(resam), function(i){var(resam[which(resam[,i]>0),i])})
  a0.imputed <- mus*(1-mus)/vars - 1
  alphas <- mean(a0.imputed)*mus
  binmat <- matrix(data=NA, nrow=N, ncol=ncol(base.psis))
  for(i in 1:ncol(binmat)){
    binmat[,i] <- sample(c(1, 0), nrow(binmat), prob=c(nonzero[i]/nrow(resam), 1 - nonzero[i]/nrow(resam)), 
                         replace=TRUE)
  }
  psi.out <- binmat
  for(i in 1:nrow(psi.out)){
    alphaset <- alphas[c(which(psi.out[i,]==1))]
    psi.out[i,c(which(psi.out[i,]==1))] <- rdirichlet(1, alphaset)
  }
  psi.out
}

mutate.grid <- function(mat){
  tochange.row <- sample(1:nrow(mat), 1)
  tochange.col <- sample(1:ncol(mat), 1)
  current <- mat[tochange.row,tochange.col]
  pool <- c(0:5)[which(!0:5 %in% current)]
  newval <- sample(pool, 1)
  changed <- mat
  changed[tochange.row,tochange.col] <- newval
  changed
}

mutate.grid2 <- function(mat, n){
  tochange.row <- sample(1:nrow(mat), n, replace=TRUE)
  tochange.col <- sample(1:ncol(mat), n)
  currents <- sapply(1:n, function(i){mat[tochange.row[i],tochange.col[i]]})
  pools <- t(sapply(1:n, function(i){c(0:5)[which(!0:5 %in% currents[i])]}))
  newvals <- t(sapply(1:n, function(i){sample(pools[i,], 1)}))
  changed <- mat
  for(i in 1:n){
    changed[tochange.row[i],tochange.col[i]] <- newvals[i]
  }
  unname(changed)
}

mutate.grid3 <- function(mat, n){
  rows <- nrow(mat)
  cols <- ncol(mat)
  asvec <- as.vector(mat)
  sampled <- sample(1:length(asvec), n, replace=FALSE)
  current <- asvec[c(sampled)]
  other <- t(sapply(1:length(current), function(i){which(!0:5 %in% current[i])}))
  template <- asvec
  for(i in 1:length(current)){
    template[i] <- sample(other[i,], 1)
  }
  matrix(template, nrow=rows, ncol=cols)
}

mutate.grid4 <- function(mat, i, dir){
  rows <- nrow(mat)
  cols <- ncol(mat)
  asvec <- as.vector(mat)
  sampled <- sample(1:length(asvec), n, replace=FALSE)
  current <- asvec[c(sampled)]
  other <- t(sapply(1:length(current), function(i){which(!0:5 %in% current[i])}))
  template <- asvec
  for(i in 1:length(current)){
    template[i] <- sample(other[i,], 1)
  }
  matrix(template, nrow=rows, ncol=cols)
}

scramble <- function(mat, psi, z){
  
}

mutate.psi <- function(psivec, delta){
  psi <- psivec
  pick <- sample(1:length(psi), 1)
  other <- which(!1:length(psi) %in% pick)
  if(psi[pick]==0){
    sign <- 1
  }
  if(psi[pick]==1){
    sign <- -1
  }
  if(psi[pick] + delta < 0){
    delta <- psi[pick]
  }
  if(psi[pick] + delta > 1){
    delta <- 1 - psi[pick]
  }
  psi[pick] <- psi[pick] + delta
  adjustment <- (-1)*delta/(length(psi) - 1)
  psi[other] <- psi[other] + adjustment
  psi[which(psi<0)] <- 0
  psi[which(psi>1)] <- 1
  psi
}

mutate.psi2 <- function(psi, index, delta, dir){
  psivec <- psi
  if(psivec[index] + dir*delta < 0){
    change <- psivec[index]*dir
  }else if(psivec[index] + dir*delta > 1){
    change <- 1 - psivec[index]
  }else{
    change <- delta*dir
  }
  psivec[index] <- psivec[index] + change
  for(i in which(1:length(psi) %in% index)){
    psivec[i] <- psivec[i] + -1*dir*change/(length(psivec)-1)
  }
  psivec
}

###Assessment:
findBest.mat <- function(matpop, truemat){
  diffs <- sapply(1:dim(matpop)[1],function(i){
    truemat <- t(rbind(truemat[,c(1:(ncol(truemat)/2))], truemat[,(ncol(truemat)/2+1):ncol(truemat)]))
    short <- nrow(matpop[i,,]) - nrow(truemat)
    if(short==1){
      zeros <- t(matrix(rep(0, ncol(truemat))))
      truemat <- rbind(truemat, zeros)
    }else if(short > 1){
      zeros <- matrix(data=0, nrow=short, ncol=ncol(truemat))
      truemat <- rbind(truemat, zeros)
    }
    sum((matpop[i,,] - truemat)^2)                                      
  })
  list(matpop[which.min(diffs),,], which.min(diffs), min(diffs))
}

findBest.vec <- function(){
  
}

findBest.psi <- function(psipop, truepsi){
  truepsi <- c(truepsi, rep(0, ncol(psipop) - length(truepsi)))
  diffs <- sapply(1:nrow(psipop),function(i){sum((psipop[i,] - truepsi)^2)})
  c(psipop[which.min(diffs),], which.min(diffs), min(diffs))
}

rank.fitness <- function(psipop, gridpop, truepsi, truemat){
  #true.resids <- 
  
}

compute.post <- function(mat, psivec, datavec, markers, theta, sigma0, ktheta, alpha, kmax,
                         transposed=FALSE){
  psi <- c(psivec, rep(0, kmax-length(psivec)))
  if(transposed==FALSE){
    cnmat <- cbind(t(mat[,1:(ncol(mat)/2)]), t(mat[,(ncol(mat)/2+1):ncol(mat)])) 
  }else{
    cnmat <- mat
  }
  blank <- matrix(rep(0, kmax*length(datavec)), nrow=kmax, ncol=length(datavec))
  for(i in 1:nrow(cnmat)){
    blank[i,] <- cnmat[i,]
  }
  cnmat <- blank
  cnp <- cnPrior(cnmat, theta)
  pp <- psiPrior(psi, alpha)
  k <- length(which(psi > 0))
  kp <- kPrior(k, ktheta)
  prior <- cnp + pp + kp
  resids <- (colSums(cnmat*psi) - datavec)^2
  ll <- sum(logLike(resids, markers, sigma0))
  prior + ll
}


###Computing accuracy
abs.acc <- function(finalmat, truemat){
  binmat <- finalmat==truemat
  length(which(as.vector(binmat)))/(nrow(binmat)*ncol(binmat))
}

computeRsquared <- function(finalmat, finalpsi, datavec){
  r2.mat <- (finalmat - truemat)^2
  r2.psi <- (finalpsi - c(truepsi, rep(0, kmax-length(truepsi))))^2
  c(r2.mat, r2.psi)
}

computePval <- function(finalmat, finalpsi, datavec, markers, sigma.cn, sigma.psi){
  sigma.mat <- sigma.cn
  r2.mat <- (finalmat - truemat)^2
  sigmas <- 
  r2.psi <- (finalpsi - truepsi)^2
  sum(dnorm(as.vector(r2.mat), 0, sigma.mat/sigmas, log=TRUE), dnorm(r2.psi, 0, sigma.psi, log=TRUE))
}

computeMetric <- function(mat, truemat){
  M <- mat %*% t(mat)
  TM <- truemat %*% t(truemat)
  diffmat <- M - MT
  eigenvals <- eigen(diffmat)$value
  sum(eigenvals)
}

runAlg <- function(FUN, params){
  #FUN(params)
  #do.call(FUN, params)
}

assess <- function(FUN, final, true){
  
}

smooth.grid <- function(vecs, posts, range){
  
}

smooth.vec <- function(vec, posts, range){
  
}

smooth.psi <- function(psis, posts){
  
}

