#Method of moments to reparametrize beta distribution.
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  c('alpha' = abs(alpha), 'beta' = abs(beta))
}

logLike <- function(resids, nmarks, sigma0){
  dnorm(resids, 0, sigma0/(nmarks)^.5, log=TRUE)
}
#Length of resids = length of nmarks

psiPrior <- function(psi, alpha=.5, kmax=5, minim=.001){
  psivec <- psi
  psivec[which(psivec<minim)] <- 0
  psivec <- sort(psivec/sum(psivec), decreasing=TRUE)
  if(length(which(psivec>0))==1){
    psivec[1] <- 1 - minim
    psivec[2] <- minim
  }
  p <- log(ddirichlet(psivec[c(which(psivec>0))], rep(alpha, length(which(psivec>0)))))
  p
}

psiPrior2 <- function(psi, alpha, kmax=5, minim=.001){
  psivec <- psi
  psivec[which(psivec<minim)] <- minim
  psivec[which(psivec>1-minim)] <- 1 - minim*length(which(psivec==minim))/length(which(psivec>1-minim))
  psivec <- psivec/sum(psivec)
  p <- log(ddirichlet(psivec, rep(alpha, length(psivec))))
  p
}

cnPrior <- function(cnmat, psi, theta){
  nonzero <- which(psi>0)
  if(class(cnmat)=='numeric' | class(cnmat)=='integer'){
    output <- sum(dgeom(abs(1 - as.vector(cnmat[c(nonzero)])), theta, log=TRUE))
  }else{
    output <- sum(dgeom(abs(1 - as.vector(cnmat[,c(nonzero)])), theta, log=TRUE)) 
  }
  output
}

cnPrior.vec <- function(cnvec, psi, theta){
  nonzero <- which(psi>0)
  sum(dgeom(abs(1 - cnvec[nonzero]), theta, log=TRUE))
}

#Prior on number of subclones
kPrior <- function(k, ktheta){
  dgeom(k, ktheta, log=TRUE)
}

#Prior on number of unique segments:
sPrior <- function(cnmat, stheta){
  dgeom(nrow(unique(cnmat)), prob=stheta)
}

priorGen <- function(psis, cnmodels, pars, cnmax=5){
  #cnpriors <- dgeom(0:cnmax, pars$theta, log=TRUE)
  #kpriors <- dgeom(1:ncol(psis), pars$ktheta, log=TRUE)
  mat <- matrix(NA, nrow=nrow(psis), ncol=nrow(cnmodels))
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      mat[i,j] <- log(ddirichlet(psis[i,][which(psis[i,]>pars$thresh)],rep(pars$alpha,
       length(which(psis[i,]>pars$thresh))))) + sum(dgeom(abs(1-cnmodels[j,][which(psis[i,]>pars$thresh)]),pars$theta, log=TRUE)) 
      + dgeom(length(which(psis[i,]>pars$thresh)),pars$ktheta,log=TRUE) + dgeom(length(unique(cnmodels[j,])), pars$ntheta, log=TRUE)
    }
  }
  mat
}

#Crossover operator for genetic algorithm
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
  if(class(parents1)=='matrix'){
    par1 <- parents1
    par2 <- parents2
    n1 <- sample(1:nrow(par1), 1)
    mats <- lapply(1:length(1:O), function(j){
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
  }else{
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
  }
  mats
}

findPsi <- function(datavec, AB){
  AB <- as.matrix(AB)
  M <- AB
  keep <- which(sapply(1:ncol(M), function(i){length(which(M[,i]==0)) < nrow(M)}))
  M <- as.matrix(M[,c(keep)])
  dat <- as.matrix(datavec)
  psi <- as.vector(solve(t(M) %*% M) %*% t(M) %*% dat)
  psi <- abs(psi)
  psi <- psi/sum(psi)
  psi <- c(psi, rep(0, ncol(AB) - length(psi)))
  psi
}
#(M(T)M)^(-1)*M(T)*W; W = c(X, Y) 

#Convert probabilities that R thinks are 0 into something usable (using logs)
convert <- function(log.probs, log=TRUE){
  if(log){
    probs.temp <- exp(log.probs)
  }else{
    probs.temp <- log.probs
  }
  if(length(which(probs.temp>0))==length(log.probs)){
    probs <- probs.temp
    probs <- probs/sum(probs)
  }else{
    maxim <- max(log.probs)
    probs <- log.probs - maxim
    probs <- exp(probs)
    minim <- min(na.omit(probs))
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
  probs[which(is.na(probs))] <- min(na.omit(probs))
  probs
}

filter <- function(data, threshold){
  indices <- which(abs(data$X - 1)>=threshold | abs(data$Y - 1)>=threshold)
  list('mat'=data[indices,], 'indices'=indices)
}

filter.mut <- function(mutdata, mu, threshold){
  indices <- which(abs(mutdata$refCounts - mu)>=threshold | abs(mutdata$varCounts - mu)>=threshold)
  ids <- mutdata$mut.id[indices]
  list('mat'=mutdata[indices,], 'ids'=ids)
}

mut.cluster <- function(muts, kmax, theta, sigma0){
  ind <- which(colnames(muts) %in% c('refCounts', 'varCounts'))
  mat <- as.matrix(muts[,ind])
  if(nrow(mat)>=2){
    d <- dist(mat)
    dmat <- matrix(nrow=nrow(mat), ncol=nrow(mat))
    for(i in 1:nrow(mat)){
      for(j in 1:nrow(mat)){
        dmat[i,j] <- sum((mat[i,] - mat[j,])^2)^.5
      }
    }
    hc <- hclust(d)
    stats <- sapply(1:kmax,function(x){
      tree <- cutree(hc, k=x)
      mean(dmat[which(tree==x),which(tree==x)])/mean(dmat)
    })
    priors <- dgeom(1:kmax, theta)
    priors <- priors/sum(priors)
    k <- which.min(priors*(stats/sum(stats)))
    pick <- cutree(hc, k=k)
    clusters <- unique(pick)
    df <- t(sapply(1:length(clusters),function(j){
      c(mean(mat[which(pick==clusters[[j]]),1]), 
        mean(mat[which(pick==clusters[[j]]),2]), 
        sd(mat[which(pick==clusters[[j]]),1]), 
        sd(mat[which(pick==clusters[[j]]),2]),
        length(which(pick==clusters[[j]])),
        j)
    }))
  }else if(nrow(mat)==1){
    df <- t(as.matrix(c(mat[,1], mat[,2], 1/(sigma0^.5), 1/(sigma0^.5), nrow(mat), 1)))
  }else{
    df <- t(as.matrix(rep(NA, 6)))
  }
  colnames(df) <- c('refCounts.mu','varCounts.mu','refCounts.sd','varCounts.sd','mutations','cluster')
  if(nrow(mat>1)){
    colors <- c('black', 'gray', 'red', 'blue', 'green', 'chartreuse',
                'cadetblue1', 'darkorange', 'cyan', 'darkgoldenrod1')
    plot(mat[,1], mat[,2], pch=16, cex=.8, col=colors[clusters],
         xlab='Ref Counts', ylab='Var Counts', xlim=c(0, 100), ylim=c(0, 100))
  }
  list('df'=as.data.frame(df), 'cluster.index'=pick)
}

#N must be less than or equal to nrow(psipool)
filter2 <- function(psipool, vecpool, datavec){
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

#Resampling matrices from a distributions
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
  for(i in 1:nrow(resam)){
    resam[i,] <- as.vector(rdirichlet(1, resam[i,]))
  }
  resam
}

resample.grid2 <- function(cnmat, posts, N){
  probs <- convert(posts)
  resam.grid <- gridmat[sample(1:nrow(gridmat), prob=probs, replace=TRUE),]
  mus <- colMeans(resam.grid)
  mus <- round(mus)
  #s2s <-
  from.dist <- t(sapply(1:length(mus), function(i){rnbinom(N, mus[i], size=r)}))
  from.dist
}

#Assumes error for copy number values estimates is symmetric
resample.grid3 <- function(cnmat, posts, N){
  probs <- convert(posts)
  resam.grid <- cnmat[sample(1:nrow(cnmat), N, prob=probs, replace=TRUE),]
  modes <- colMeans(resam.grid)
  mu.devs <- sapply(1:ncol(resam.grid), function(i){mean(abs(resam.grid[,i] - modes[i]))})
  devs <- sapply(1:length(mu.devs), function(i){
    abs.devs <- rgeom2(N, mu.dev=mu.devs[i])
    signs <- sample(c(-1, 1), N, replace=TRUE)
    abs.devs*signs
  })
  out <- resam.grid + devs
  asvec <- as.vector(out)
  asvec[which(asvec>5)] <- 5
  asvec[which(asvec<0)] <- 0
  matrix(asvec, nrow=nrow(out), ncol=ncol(out))
}

#Mutation operator
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
  for(i in 1:length(sampled)){
    template[sampled[i]] <- sample(c(0:5)[other[i,]], 1)
  }
  newmat <- matrix(template, nrow=rows, ncol=cols)
  newmat
}

mutate.grid4 <- function(mat, i, dir){
  rows <- nrow(mat)
  cols <- ncol(mat)
  asvec <- as.vector(mat)
  sampled <- sample(1:length(asvec), n, replace=FALSE)
  current <- asvec[c(sampled)]
  other <- t(sapply(1:length(current), function(j){which(!0:5 %in% current[j])}))
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
  sort(psi, decreasing=TRUE)
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

mutate.psi3 <- function(psivec, deltarange){
  delta <- runif(1, 0, deltarange)
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
  psi <- psi/sum(psi)
  sort(psi, decreasing=TRUE)
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

findBest.psi <- function(psipop, truepsi){
  truepsi <- c(truepsi, rep(0, ncol(psipop) - length(truepsi)))
  diffs <- sapply(1:nrow(psipop),function(i){sum((psipop[i,] - truepsi)^2)})
  c(psipop[which.min(diffs),], which.min(diffs), min(diffs))
}


compute.post <- function(mat, psivec, datavec, markers, theta, sigma0, ktheta, alpha, kmax, transposed=FALSE){
  if(length(psivec)<kmax){
    psi <- c(psivec, rep(0, kmax-length(psivec)))
  }else{
    psi <- psivec
  }
  if(transposed==FALSE){
    cnmat <- unname(t(mat))
  }else{
    cnmat <- mat
  }
  blank <- matrix(rep(0, kmax*length(datavec)), nrow=kmax, ncol=length(datavec))
  for(i in 1:nrow(cnmat)){
    blank[i,] <- cnmat[i,]
  }
  cnmat <- blank
  psi <- unlist(psi)
  psi <- psi/sum(psi)
  cnp <- cnPrior(cnmat, psi, theta)
  pp <- psiPrior(unlist(psi), alpha)
  k <- length(which(unlist(psi) > 0))
  kp <- kPrior(k, ktheta)
  prior <- cnp + pp + kp
  resids <- (colSums(cnmat*unlist(psi)) - datavec)^2
  ll <- sum(logLike(resids, markers, sigma0))
  post <- prior + ll
  post
}


###Computing accuracy metrics
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

smooth.grid <- function(vecs, posts, range){
  
}

smooth.vec <- function(vec, posts, range){
  
}

smooth.psi <- function(psis, posts){
  
}

compDist2 <- function(subj.vec, ref.vec, subj.psi, ref.psi, markers, scale=1, x0=1){
  diffmat <- ref.vec - subj.vec
  psidiff.sq <- (ref.psi - subj.psi)^2
  weighted.diff <- ref.psi*t(diffmat)
  rmat.sq <- ((markers/sum(markers))*t(weighted.diff))^2
  s.cn <- (sum((psidiff.sq + x0)*t(rmat.sq)))
  s.psi <- sum(abs(psidiff.sq))
  s <- (s.cn/scale + s.psi)^.5
  s
}

#Current working version
compDist3 <- function(subj.vec, ref.vec, subj.psi, ref.psi, markers, scale=1){
  diffmat <- ref.vec - subj.vec
  psidiff.sq <- (ref.psi - subj.psi)^2
  weighted.diff <- ref.psi*t(diffmat)
  rmat.sq <- ((markers/sum(markers))*t(weighted.diff))^2
  s.cn <- (sum(rmat.sq))
  s.psi <- sum(abs(psidiff.sq))
  s <- (s.cn/scale + s.psi)^.5
  s
}

compEditDist <- function(truth, imputed, grain, criteria='absolute'){
  if(criteria=='absolute'){
    
  }else{
    
  }
}

success <- function(mat, truemat, psi, truepsi, psi.con, mat.con){
  length(which(as.vector(mat==truemat)))/(nrow(mat)*ncol(mat))
}

#Optionally weighted
compMatDist <- function(mat, truemat, markers=NULL){
  diffmat <- abs(mat - truemat)
  if(is.null(markers)){
    weights <- rep(1, nrow(mat))
  }else{
    weights <- markers/sum(markers)
  }
  weighted <- weights*diffmat
  (sum(weighted^2))^.5
}

compPsiDist <- function(psi, truepsi){
  sum((sort(psi) - sort(truepsi))^2)^.5
}

findNeighbors <- function(mat, dist, markers, inclusive=TRUE, mindist=0, max.cn=5){
  wts <- markers/sum(markers)
  dist.temp <- min(dist, sum(as.vector(sapply(1:length(wts),function(i){wts[i]*(rep(1,ncol(mat)))}))^2)^.5)
  if(inclusive){
    distance <- runif(1, mindist, dist.temp)
  }else{
    distance <- dist.temp
  }
  set <- rep(NA, nrow(mat)*ncol(mat))
  val <- 0
  iter <- 1
  while(val < distance) {
    pool <- which(sapply(1:nrow(mat), function(i){length(which(set==i))<ncol(mat)}))
    oldset <- set
    if(length(pool)>1){
      set[iter] <- sample(pool, 1)
    }else{
      set[iter] <- pool 
    }
    oldval <- val
    inc <- sapply(1:nrow(mat), function(i){length(which(set==i))})
    iter <- iter + 1
    val <- sum(wts[set[which(!is.na(set))]]^2)^.5
  }
  if(abs(oldval - distance) < abs(val - distance)){
    iter <- iter - 1
  }
  set <- set[1:iter]
  rowN <- sapply(1:nrow(mat), function(i){length(which(set==i))})
  deltamat <- matrix(0, nrow=nrow(mat), ncol=ncol(mat))
  for(i in which(rowN>0)){
    cols <- sample(1:ncol(mat), rowN[i], replace=FALSE)
    for(j in cols){
      if(mat[i,j]==max.cn){
        deltamat[i,j] <- -1
      }else if(mat[i,j]==0){
        deltamat[i,j] <- 1
      }else{
        deltamat[i,j] <- sample(c(-1, 1), 1)
      }
    }
  }
  mat + deltamat
}

findNeighbors2 <- function(mat, D, markers, inclusive=TRUE, max.cn=5){
  mindist <- 1/nrow(mat)
  if(inclusive){
    dist <- runif(1, mindist, D)
  }else{
    dist <- D
  }
  temp <- dist^2*(nrow(mat)^2)
  asvec <- as.vector(mat)
  x <- round(temp)
  x <- min(x, length(asvec))
  dist <- (x*(1/(nrow(mat)))^2)^.5
  indices <- sample(1:length(asvec), x, replace=FALSE)
  for(i in 1:x){
    curr <- asvec[indices[i]]
    if(curr==max.cn){
      asvec[indices[i]] <- curr - 1
    }else if(curr==0){
      asvec[indices[i]] <- curr + 1
    }else{
      asvec[indices[i]] <- curr + sample(c(1, -1), 1)
    }
  }
  matrix(asvec, nrow=nrow(mat), ncol=ncol(mat))
}

findNeighbors.psi <- function(psi, dist){
  
}

#Compute direction of change from one parameter set to another
compDir <- function(subj.vec, ref.vec, subj.psi, ref.psi){
  diffmat <- ref.vec - subj.vec
  psidiff <- ref.psi - subj.psi
  weighted.diff <- psidiff*t(diffmat)
  as.vector(weighted.diff)
}

compAngle <- function(start.psi, start.cn, changed.psi, changed.cn, true.psi, true.cn){
  diff.cn1 <- true.cn - start.cn
  diff.psi1 <- true.psi - start.psi
  dir1 <- as.vector(diff.cn1*diff.psi1)
  diff.cn2 <- true.cn - changed.cn
  diff.psi2 <- true.psi - changed.psi
  dir2 <- as.vector(diff.cn2*diff.psi2)
  dot.prod <- sum(dir1*dir2)
  norm.prod <- (sum(dir1^2)*sum(dir2^2))^.5
  cos.a <- dot.prod/norm.prod
  180*acos(cos.a)/pi
}

compTrace <- function(arguments){
  res <- arguments[[1]]
  true <- arguments[[2]]
  final <- res[[1]]
  mat <- final[[1]]
  psi <- final[[2]]
  truemat <- true[[1]]
  truemat <- rbind(truemat[,1:(ncol(truemat)/2)], truemat[,(ncol(truemat)/2+1):ncol(truemat)])
  blank <- matrix(0, nrow=nrow(mat), ncol=ncol(mat))
  blank[,c(1:nrow(truemat))] <- truemat
  truemat <- blank
  M <- mat %*% t(mat)
  MT <- truemat %*% t(truemat)
  diffmat <- M - MT
  eigenvals <- eigen(diffmat)$value
  sum(eigenvals)
}

oldRunAlg <- function(FUN, params){
  #FUN(params)
  do.call(FUN, params)
}

assess <- function(FUN, arguments){
  FUN(arguments)
}

##Some probabity distribution functions:
rgeom2 <- function(N, mu.dev){
  p <- 1/(mu.dev + 1)
  rgeom(N, prob=p)
}

dgeom2 <- function(x, mu.dev, log=FALSE){
  p <- 1/(mu.dev + 1)
  dgeom(x, prob=p, log=log)
}

rbeta2 <- function(n, mu, sigma){
  params <- estBetaParams(mu, sigma^2)
  rbeta(n, params[1], params[2])
}

dbeta2 <- function(x, mu, sigma, log=FALSE){
  params <- estBetaParams(mu, sigma^2)
  dbeta(x, params[1], params[2], log=log)
}

#scaled beta
rsb <- function(n, alpha, beta, scalar, minim){
  betas <- rbeta(n, alpha, beta)
  scaled <- round(betas*scalar)
  scaled + minim
}

dsb <- function(x, alpha, beta, scalar, minim, log=FALSE){
  betas <- (x - minim)/scalar
  dbeta(betas, alpha, beta, log=log)
}

rgamma2 <- function(n, mu, sigma){
  k <- mu^2/sigma^2
  theta <- sigma^2/mu
  rgamma(n, shape=k, scale=theta)
}

dgamma2 <- function(x, mu, sigma, log=FALSE){
  k <- mu^2/sigma^2
  theta <- sigma^2/mu
  dgamma(x, shape=k, scale=theta, log=log)
}

#Discrete uniform sum cumulative distribution


##Test to see if true parameter set has best posterior probability
testpost <- function(pars, n, q, r, bestpost, mainlab, truemat, truepsi){
  #q: max number of mutations; at most, nrow*ncol; q <= length(colorset)
  #colorset <- c('cadetblue1', 'chartreuse', 'chocolate1', 'coral3', 'cornflowerblue', 'cyan', 'blue',
  #              'brown1', 'darkcyan', 'darkgoldenrod2', 'darkgreen', 'darkmagenta')
  datavec <- pars$datavec
  markers <- pars$markers
  tpost <- compute.post(truemat, truepsi, datavec=datavec, markers=markers, 
                        theta=pars$theta, sigma0=pars$sigma0, ktheta=pars$ktheta, alpha=pars$alpha, kmax=pars$kmax,
                        transposed=FALSE)
  jset <- sapply(1:q, function(j){rep(j, n)})
  mset <- lapply(1:(q*n), function(j){
    mutate.grid3(truemat, jset[j])
  })
  mset <- rep(mset, r+1)
  deltas <- runif(r, min=0, max=.2)
  mpsiset <- lapply(1:(r+1), function(i){
    if(i==1){
      out <- truepsi
    }else{
      out <- mutate.psi(truepsi, deltas[i-1])
    }
    out
  })
  mpsiset <- lapply(1:length(mpsiset), function(i){rep(mpsiset[i], q*n)})
  mpsiset <- unlist(mpsiset, recursive <- FALSE)
  mposts <- sapply(1:length(mset), function(j){
    mat <- mset[[j]]
    if(nrow(mat)==length(datavec)/2){
      mat <- rbind(mat[,1:(ncol(mat)/2)], mat[,(ncol(mat)/2+1):ncol(mat)])
    }
    compute.post(mat=mat, psivec=mpsiset[[j]], datavec=datavec, markers=markers, 
                 theta=pars$theta, sigma0=pars$sigma0, ktheta=pars$ktheta, alpha=pars$alpha, kmax=pars$kmax)
  })
  dists <- sapply(1:length(mpsiset), function(i){
    mat <- mset[[i]]
    if(nrow(mat)==length(datavec)/2){
      mat <- rbind(mat[,1:(ncol(mat)/2)], mat[,(ncol(mat)/2+1):ncol(mat)])
    }
    mat <- cbind(mat, matrix(0, nrow=nrow(mat), ncol=pars$kmax - ncol(mat)))
    compDist3(mat, truemat, mpsiset[[i]], truepsi, markers)
  })
  #colors <- sapply(1:q,function(i){rep(colorset[i], n)})
  #colors <- c('black', rep(colors, r+1))
  colfunc <- colorRampPalette(c('hotpink', 'purple'))
  #plot(rep(1,5),col=colfunc(5),pch=19,cex=3)
  colors <- colfunc(5)
  colors <- rev(colors)
  borders <- (0:4)/5
  colors <- sapply(dists, function(i){colors[max(which(borders<=i))]})
  colors <- c('black', colors)
  y <- c(tpost, mposts)
  x <- 1:length(y)
  sizes <- rep(1, length(x))
  sizes[1] <- 2
  #if(max(y)==tpost){
  #  mainline <- 'Max post. = true params'
  #}else{
  #  mainline <- 'Max post. =/= true params'
  #}
  minim <- min(c(y, bestpost))
  maxim <- max(c(y, tpost))
  plot(x, y, ylim=c(minim, maxim), col=colors, pch=16, cex=sizes, main=mainlab)
  #'posts'=y
}

#Evaluate algorithm results
eval1 <- function(res, alg, pars, alone=TRUE, truepsi, truemat){
  kmax <- pars$kmax
  markers <- pars$markers
  datavec <- pars$datavec
  final <- res[[1]]
  arrays <- res[[2]]
  cnarray <- arrays[[1]]
  psiarray <- arrays[[2]]
  postmat <- arrays[[3]]
  llarray <- arrays[[4]]
  exparray <- arrays[[5]]
  iters <- pars$iters
  #true.psi <- truepsi
  #true.psi <- c(true.psi, rep(0, kmax - length(true.psi)))
  #tmat <- rbind(truemat[,(1:(ncol(truemat)/2))], truemat[,(ncol(truemat)/2+1):ncol(truemat)])
  #true.mat <- tmat
  #true.mat <- cbind(true.mat, matrix(data=0, ncol=kmax - ncol(true.mat), nrow=nrow(true.mat)))
  #order.psi <- sort(true.psi, decreasing=TRUE, index.return=TRUE)$ix
  #true.psi <- true.psi[order.psi]
  #true.mat <- true.mat[,c(order.psi)]
  tpost <- compute.post(truemat, truepsi, datavec, markers, theta=pars$theta, sigma0=pars$sigma0, 
                        ktheta=pars$ktheta, alpha=pars$alpha, kmax=5, transposed=FALSE)
  tll <- sum(logLike((colSums(unlist(truepsi)*t(truemat)) - datavec)^2, markers, pars$sigma0))
  measurements <- as.data.frame(t(sapply(1:iters, function(i){
    if(alg=='mcmc'){
      min <- (i-1)*pars$N + 1
      max <- i*pars$N
      premax <- min-1
      premin <- min - pars$N
      if(i>1){
        posts.old <- rowSums(postmat[premin:premax,])
        oldbestpsi <- psiarray[,premin:premax][,which.max(posts.old)]
        oldbestmat <- cnarray[,,premin:premax][,,which.max(posts.old)]
      }
      posts.new <- rowSums(postmat[min:max,])
      maxpostnew <- max(posts.new)
      maxllnew <- max(rowSums(llarray[min:max,]))
      newbestpsi <- psiarray[,min:max][,which.max(posts.new)]
      newbestmat <- cnarray[,,min:max][,,which.max(posts.new)]
      exps <- exparray[,min:max][,which.max(posts.new)]
    }else if(alg=='ga'){
      imax <- which.max(postmat[,i])
      posts.new <- postmat[,i]
      maxpostnew <- max(posts.new)
      maxllnew <- max(llarray[,i])
      newbestmat <- cnarray[imax,,,i]
      newbestpsi <- as.vector(psiarray[imax,,i])
      if(i>1){
        posts.old <- postmat[,i-1]
        imax.minus <- which.max(posts.old)
        oldbestpsi <- as.vector(psiarray[imax.minus,,i-1])
        oldbestmat <- cnarray[imax.minus,,,i-1]
      }
      exps <- exparray[imax,,i]
    }else{
      maxpostnew <- postmat[i]
      maxllnew <- llarray[i]
      newbestmat <- t(cnarray[i,,])
      newbestpsi <- as.vector(psiarray[i,])
      if(i>1){
        oldbestmat <- t(cnarray[i-1,,])
        oldbestpsi <- as.vector(psiarray[i-1,])
      }
      exps <- exparray[i,]
    }
    if(i>1){
      ang <- compAngle(oldbestpsi, oldbestmat, newbestpsi, newbestmat, truepsi, truemat)
      delta <- compDist3(newbestmat, oldbestmat, newbestpsi, oldbestpsi, markers)
    }else{
      ang <- NA
      delta <- NA
    }
    dist <- compDist3(newbestmat, truemat, newbestpsi, truepsi, markers)
    postpick <- maxpostnew
    llpick <- maxllnew
    exp.resid <- sum(abs(exps - datavec))
    #trace <- compTrace()
    c(ang, dist, delta, postpick, llpick, exp.resid)
  })))
  colnames(measurements) <- c('Angle', 'Distance', 'Delta', 'Post', 'll', 'exp.resids')
  if(alone==TRUE){
    par(mfrow=c(2,3))
  }
  
  #plot(2:nrow(measurements), measurements$Angle[2:nrow(measurements)], xlim=c(1, nrow(measurements)),
  #     type='l', xlab='iteration', ylab='Angle')
  #abline(h=c(0, 60, 180), col='black')
  plot(measurements$ll, ylim=c(min(c(measurements$ll, tll)), max(c(measurements$ll, tll))), 
       type='l', xlab='iteration', ylab='Log Likelihood')
  abline(h=tll, col='red')
  #plot(measurements$exp.resids/length(datavec), type='l', ylim=c(min(c(measurements$exp.resids, 0)), 
  #    max(measurements$exp.resids)), xlab='iteration', ylab='Exp. Resid.')
  true.exp <- colSums(sort(truepsi, decreasing=TRUE)*t(truemat))
  #abline(h=sum((true.exp-datavec)^2)/length(datavec), col='red')
  plot(measurements$Distance, type='l', xlab='iteration', ylab='Distance')
  abline(h=0, col='red')
  #plot(2:nrow(measurements), measurements$Delta[2:nrow(measurements)], xlim=c(1, nrow(measurements)),
  #     type='l', xlab='iteration', ylab='Delta')
  #abline(h=, )
  plot(measurements$Post, ylim=c(min(c(measurements$Post, tpost)), max(c(measurements$Post, tpost))), 
       type='l', xlab='iteration', ylab='Post')
  abline(h=tpost, col='red')
  bestpost <- max(measurements$Post)
  if(bestpost > tpost){
    mainlab <- 'Max post. > true post.'
  }else{
    mainlab <- 'Max post. = true post.'
  }
  testpost(pars, 8, 10, 12, bestpost, mainlab=mainlab, truemat=truemat, truepsi=truepsi)
  abline(h=bestpost, col='black')
  abline(h=tpost, col='red')
  plotfun(datavec, colSums(final[[2]]*t(final[[1]])), markers, alone=FALSE)
  measurements
}

eval2 <- function(truth, data, res, pars, tcn, thresh, z=FALSE){
  if(z){
    Euclidean <- NA
    meanResidSquared <- NA
    psidist <- compPsiDist(psi=res$psi, truepsi=truth$truepsi)
    cndist <- NA
    psiDiff <- sort(res$psi) - sort(truth$truepsi)
    cnDiff <- NA
    cnDiff.avg <- NA
    cnDiff.var <- NA
    #acc.rate <- sum(sapply(1:length(truth$truepsi), function(k){length(which(truth$truearray[,k,]==res$array[,,k]))}))/(length(as.vector(truth$truearray)))
    acc.rate <- NA
    cn.total <- matrix(NA, nrow=dim(res$array)[1], ncol=dim(res$array)[3])
    for(i in 1:ncol(cn.total)){
      for(j in 1:nrow(cn.total)){
        cn.total[j,i] <- tcn[which(res$array[j,,i]==1)]
      }
    }
    cntrue.total <- matrix(NA, nrow=dim(truth$truearray)[1], ncol=dim(truth$truearray)[2])
    for(i in 1:ncol(cntrue.total)){
      for(j in 1:nrow(cntrue.total)){
        cntrue.total[j,i] <- tcn[which(truth$truearray[j,i,]==1)]
      }
    }
    tot.cndiff <- as.vector(cn.total - cntrue.total)
    tot.cndiff.avg <- mean(tot.cndiff)
    tot.cndiff.var <- var(tot.cndiff)
    tot.acc <- length(which(as.vector(cntrue.total==cn.total)))/(nrow(cn.total)*ncol(cn.total))
    kdiff <- length(which(res$psi>thresh)) - length(which(truth$truepsi>0))
    runtime <- NA
  }else{
    Euclidean <- compDist3(subj.vec=res$mat, ref.vec=truth$truemat, subj.psi=res$psi, ref.psi=truth$truepsi, markers=data$markers)
    meanResidSquared <- mean((data$datavec - colSums(res$psi*t(res$mat)))^2)^.5
    psidist <- compPsiDist(psi=res$psi, truepsi=truth$truepsi)
    cndist <- compMatDist(mat=res$mat, truemat=truth$truemat)
    psiDiff <- sort(res$psi) - sort(truth$truepsi)
    cnDiff <- as.vector(res$mat[,sort(res$psi, index.return=TRUE)$ix] - truth$truemat[,sort(truth$truepsi, index.return=TRUE)$ix])
    cnDiff.avg <- mean(cnDiff)
    cnDiff.var <- var(cnDiff)
    cn.total <- res$mat[1:(nrow(res$mat)/2),] + res$mat[(nrow(res$mat)/2+1):nrow(res$mat),]
    acc.rate <- length(which(as.vector(truth$truemat==res$mat)))/(length(as.vector(truth$truemat)))
    tot.cndiff <- as.vector(truth$truemat[1:(nrow(truth$truemat)/2),] + truth$truemat[(nrow(truth$truemat)/2+1):nrow(truth$truemat),] - 
      cn.total)
    tot.cndiff.avg <- mean(tot.cndiff)
    tot.cndiff.var <- var(tot.cndiff)
    tot.acc <- length(which(as.vector(truth$truemat[1:(nrow(truth$truemat)/2),] + 
      truth$truemat[(nrow(truth$truemat)/2+1):nrow(truth$truemat),] == cn.total)))/(nrow(cn.total)*ncol(cn.total))
    kdiff <- length(which(res$psi>thresh)) - length(which(truth$truepsi>0))
    runtime <- res$runtime
    if(runtime<0){
      runtime <- -runtime
    }
  }
  list('Euclidean'=Euclidean,'meanResidSquared'=meanResidSquared,'psidist'=psidist,'cndist'=cndist,'cnDiff.avg'=cnDiff.avg,
       'cnDiff.var'=cnDiff.var,'acc.rate'=acc.rate, 'psidiff1'=psiDiff[1], 'psidiff2'=psiDiff[2], 'psidiff3'=psiDiff[3],
       'psidiff4'=psiDiff[4], 'psidiff5'=psiDiff[5], 'tot.cndiff.avg'=tot.cndiff.avg, 'tot.cndiff.var'=tot.cndiff.var, 
       'tot.acc'=tot.acc, 'kdiff'=kdiff, 'runtime'=runtime)
}

visualize <- function(res){
  
}

newpopGen <- function(cn, psis, posts, n, include){
  N <- dim(cn)[1]
  nset <- round(runif(n, 1, 10))
  picks <- sample(1:length(posts), n, replace=TRUE, prob=convert(posts))
  mats <- lapply(picks, function(j){cn[j,,]})
  vecs <- t(sapply(picks, function(j){psis[j,]}))
  mats <- lapply(1:length(mats), function(j){t(mutate.grid3(mats[[j]], nset[j]))})
  vecs <- t(sapply(1:nrow(vecs), function(j){
    newpsi <- mutate.psi3(vecs[j,], .1)
    newpsi/sum(newpsi)
  }))
  list(mats, vecs)
}

plotfun <- function(cndat, cn.imputed, markers, alone=TRUE, mainlab=NULL){
  Xdat <- cndat[1:(length(cndat)/2)]
  Ydat <- cndat[(length(cndat)/2+1):(length(cndat))]
  X.imputed <- cn.imputed[1:(length(cndat)/2)]
  Y.imputed <- cn.imputed[(length(cndat)/2+1):(length(cndat))]
  inner <- function(data, color){
    coords <- lapply(1:length(data),function(i){
      if(i==1){
        xstart <- 1
        xend <- markers[1]
      }else{
        xstart <- (sum(markers[1:(i-1)])+1)
        xend <- sum(markers[1:i]) 
      }
      y <- data[i]
      #segments(x0=xstart, y0=y, x1=xend, y1=y, col=color)
      x <- c(xstart, xend)
      y <- rep(y, 2)
      list(x, y)
    })
    coords <- data.frame('x'=unlist(lapply(1:length(coords), function(i){coords[[i]][[1]]})),
                         'y'=unlist(lapply(1:length(coords), function(i){coords[[i]][[2]]})))
    lines(coords[,1], coords[,2], col=color)
  }
  if(alone){
    par(mfrow=c(2,1))
  }
  if(!is.null(mainlab)){
    main1 <- mainlab
  }else{
    main1 <- "Data vs. Imputed CN (X)"
  }
  plot(1, type="n", xlab="Position", ylab="Imputed Copy Number (X)", xlim=c(0, sum(markers[1:(length(cndat)/2)])), 
       ylim=c(0, 4), main=main1)
  inner(Xdat, color='red')
  inner(X.imputed, color='black')
  plot(1, type="n", xlab="Position", ylab="Imputed Copy Number (Y)", xlim=c(0, sum(markers[1:(length(cndat)/2)])), 
       ylim=c(0, 4), main="Data vs. Imputed CN (Y)")
  inner(Ydat, color='red')
  inner(Y.imputed, color='black')
}

cnplotfun <- function(data, chr, mainlab){
  dat <- data[data$chr==chr,]
  ylabel <- paste('Copy Number (', chr, ')', sep='')
  par(mfrow=c(2,1))
  plot(0, xlim=c(0, max(dat$end)), ylim=c(0, 2), main=mainlab, xlab='Position', ylab=ylabel)
  segments(x0=dat$start, x1=dat$end, y0=dat$X, y1=dat$X)
  plot(0, xlim=c(0, max(dat$end)), ylim=c(0, 2), main=mainlab, xlab='Position', ylab=ylabel)
  segments(x0=dat$start, x1=dat$end, y0=dat$Y, y1=dat$Y)
}

mixplotfun <- function(mixturedata, data1, data2, chr){
  dat1 <- data1[data1$chr==chr,]
  dat2 <- data2[data2$chr==chr,]
  mix.dat <- mixturedata[mixturedata$chr==chr,]
  ylabel <- paste('Copy Number (', chr, ')', sep='')
  par(mfrow=c(2,1))
  plot(0, xlim=c(0, max(dat$end)), ylim=c(0, 2), main=paste('Chr ', chr, sep=''), xlab='Position', ylab=ylabel)
  segments(x0=dat1$start, x1=dat1$end, y0=dat1$X, y1=dat1$X, col='red')
  segments(x0=dat2$start, x1=dat2$end, y0=dat2$X, y1=dat2$X, col='blue')
  segments(x0=mix.dat$start, x1=mix.dat$end, y0=mix.dat$X, y1=mix.dat$X, col='black')
  plot(0, xlim=c(0, max(dat$end)), ylim=c(0, 2), xlab='Position', ylab=ylabel)
  segments(x0=dat1$start, x1=dat1$end, y0=dat1$Y, y1=dat1$Y, col='red')
  segments(x0=dat2$start, x1=dat2$end, y0=dat2$Y, y1=dat2$Y, col='blue')
  segments(x0=mix.dat$start, x1=mix.dat$end, y0=mix.dat$Y, y1=mix.dat$Y, col='black')
}

compPost <- function(mat, psi, datavec, nmarks, pars){
  sigma0 <- pars$sigma0
  alpha <- pars$alpha
  theta <- pars$theta
  ktheta <- pars$ktheta
  kmax <- pars$kmax
  k <- length(which(psi>0))
  ll <- sum(logLike(resids=(colSums(psi*t(mat)) - datavec), nmarks, sigma0))
  p <- ll + cnPrior(mat, psi, theta) 
  + psiPrior(psi, alpha, kmax) + kPrior(k, ktheta)
  c('post'=p, 'll'=ll)
}

compPost2 <- function(ary, psi, data, nmarks, pars){
  k <- length(which(psi>0))
  compmat <- matrix(NA, nrow=dim(ary)[1], ncol=dim(ary)[3])
  for(i in 1:nrow(compmat)){
    for(j in 1:ncol(compmat)){
      compmat[i,j] <- which(ary[i,,j]==1)
    }
  }
  lrrs <- pars$lrrs
  bafs <- pars$bafs
  sigma0.lrr <- pars$sigma0.lrr
  sigma0.baf <- pars$sigma0.baf
  theta <- pars$theta
  ktheta <- pars$ktheta
  alpha <- pars$alpha
  cnpriors <- pars$cnpriors
  compset <- as.vector(compmat)
  ll <- sum(logLike2(psi=psi, zarray=ary, data=data, sigma0.lrr=sigma0.lrr, sigma0.baf=sigma0.baf,
                       lrrs=lrrs, bafs=bafs))
  p <- psiPrior(psi, alpha, kmax) + kPrior(k, ktheta) + sum(log(cnpriors[compset]))
  c('post'=p, 'll'=ll)
}

zero <- function(post, mat, psi, datavec, markers, pars){
  altpsi <- 
  altmat <- 
  altpost <- compPost(altmat, altpsi, datavec, markers, pars)
  if(altpost > post){
    output <- list('mat'=altmat, 'psi'=altpsi, 'post'=altpost)
  }
  else{
    output <- list('mat'=mat, 'psi'=psi, 'post'=post)
  }
  outpu
}

###Setting classes:
setClass('result', representation(finalPsi='array', finalMat='array', finalPost='array',
      psiArray='array', cnarray='array', postarray='array', params='list'))
setClass('assessment', representation(euclidean='numeric', trace='numeric'))

visualize <- function(res){
  
}

ptgen <- function(truepsi, truemat, min.dist, max.dist){
  
}

###Model class 2: Convert to compartmental model
#comps: normal: 1; LOH: 2; loss: 3; double loss: 4; gain of 1: 5; gain of 2: 6
convertSim <- function(cnmat, max.cn, comps=6){
  totalmat <- cnmat[1:(nrow(cnmat)/2),] + cnmat[(nrow(cnmat)/2 + 1):nrow(cnmat),]
  afmat <- cnmat[(nrow(cnmat)/2 + 1):nrow(cnmat),]/totalmat
  zarray <- array(0, dim=c(nrow(totalmat), ncol(totalmat), comps))
  for(i in 1:nrow(totalmat)){
    for(j in 1:ncol(totalmat)){
      cn <- totalmat[i,j]
      af <- afmat[i,j]
      if(cn==2 & af==.5){
        comp <- 1
      }else if(cn==2 & af!=.5){
        comp <- 2
      }else if(cn==1){
        comp <- 3
      }else if(cn==0){
        comp <- 4
      }else{
        comp <- min(cn + 2, 6)
      }
      zarray[i,j,comp] <- 1
    }
  }
  zarray
}

#Likely function for compartments
logLike2 <- function(psi, zarray, data, sigma0.lrr, sigma0.baf, lrrs, bafs){
  lrrmat <- sapply(1:dim(zarray)[3],function(i){sapply(1:dim(zarray)[1], 
     function(j){zarray[j,,i]*lrrs[which(zarray[j,,i]>0)]})})
  bafmat <- sapply(1:dim(zarray)[3],function(i){sapply(1:dim(zarray)[1], 
    function(j){zarray[j,,i]*bafs[which(zarray[j,,i]>0)]})})
  lrr.imputed <- colSums(psi*t(lrrmat))
  baf.imputed <- colSums(psi*t(bafmat))
  lls <- sapply(1:nrow(data), function(i){dnorm(data$lrr[i], lrr.imputed[i], sigma0.lrr/(data$markers[i])^.5, 
    log=TRUE) + dnorm(data$baf[i], baf.imputed[i], sigma0.baf/(data$markers[i])^.5, log=TRUE)})
  lls
}

editZ <- function(zarray, N, probs){
  n <- round(runif(1, 1, min(N, dim(zarray)[1])))
  clonepicks <- sample(1:(dim(zarray)[3]), n, replace=TRUE)
  rowpicks <- sample(1:(dim(zarray)[1]), n, replace=FALSE)
  newarray <- zarray
  for(i in 1:n){
    current <- which(zarray[rowpicks[i],,clonepicks[i]]==1)
    pool <- which(!1:length(probs) %in% current)
    alt <- sample(pool, 1, prob=probs[pool]/sum(probs[pool]))
    newarray[rowpicks[i],current,clonepicks[i]] <- 0
    newarray[rowpicks[i],alt,clonepicks[i]] <- 1
  }
  newarray
}

zPrior <- function(zarray, probs){
  probs <- probs/sum(probs)
  zpriors <- sapply(1:(dim(zarray)[1]), function(i){
    sum(sapply(1:dim(zarray)[3], function(j){log(probs[which(zarray[i,,j]==1)])}))
  })
  zpriors
}

compNormDist <- function(truepars, inferred, markers){
  wts <- markers/sum(markers)
  truemat <- truepars$truemat
  mat <- inferred$mat
  truepsi <- truepars$truepsi
  psi <- inferred$psi
  dist <- sum(c((psi - truepsi)^2, (as.vector(wts*t(mat - truemat))/max.cn)^2))^.5
  dist
}

pfun <- function(truepars, inferred, markers){
  wts <- markers/sum(markers)
  truemat <- truepars$truemat
  mat <- inferred$mat
  truepsi <- truepars$truepsi
  psi <- inferred$psi
}

genSimplex <- function(N, k, reps=1){
  simplex <- t(xsimplex(k, N))
  for(i in 1:nrow(simplex)){
    simplex[i,] <- sort(simplex[i,], decreasing=TRUE)
  }
  simplex <- unique(simplex)
  psis <- lapply(rep(1:nrow(simplex), reps), function(i){simplex[i,]/N})
  psis <- Reduce(rbind, psis)
  unname(psis)
}


#Resampling:
resample <- function(parsets, probs, pars, select, generate, subsample, subsample.psi,
   subsample.pool, psigen, data, D, deltarange, markers, mrate=.75, p, replace=TRUE){
  X <- subsample
  if(select=='sample'){
    templates <- parsets[sample(1:length(probs), X, prob=convert(probs), replace=replace)]
  }else if(select=='bestn'){
    templates <- parsets[sort(convert(probs), decreasing=TRUE, index.return=TRUE)$ix[1:X]]
  }
  per <- length(probs)/subsample
  temp.psis <- t(sapply(1:length(templates),function(i){templates[[i]]$psi}))
  nonzeros <- which(sapply(1:nrow(temp.psis), function(i){length(which(temp.psis[i,]==0))==0}))
  if(length(nonzeros)>=1){
    temp.psis <- temp.psis[nonzeros,]
  }else{
    temp.psis <- rdirichlet(nrow(temp.psis), alpha=c(.99, .01/(rep(ncol(temp.psis)-1, ncol(temp.psis)-1))))
  }
  if(class(temp.psis)=='numeric'){
    temp.psis <- t(matrix(rep(temp.psis, 2), nrow=2, ncol=2))
    temp.psis[2,] <- as.vector(rdirichlet(1, alpha=temp.psis[1,]))
  }
  alpha <- try(abs(dirichlet.mle(temp.psis, TRUE)$alpha), silent=TRUE)
  if(class(alpha)=='try-error'){
    alpha <- colMeans(temp.psis)
  }
  if(length(which(alpha<.01))>0){
    alpha[which(alpha<.01)] <- .01
  }
  out <- lapply(1:length(probs),function(i){
    ind <- floor((i-1)/per)+1
    base.psi <- templates[[ind]]$psi
    base.mat <- templates[[ind]]$mat
    curr.post <- compPost(base.mat, base.psi, datavec, markers, pars)[1]
    if(generate=='crossover'){
      parent1 <- templates[[ind]]$mat
      parent2 <- templates[[max(ind-1, X)]]$mat
      crossed <- crossover2(parents1=parent1, parents2=parent2, O=1)
      mat <- crossed[[1]]
      mut <- sample(c(TRUE, FALSE), 1, prob=c(mrate, 1-mrate))
      if(mut){
        n <- runif(1, 1, round(nrow(mat)*ncol(mat))/7)
        mat <- mutate.grid3(mat, n)
      }
    }else if(generate=='mutate'){
      mat <- templates[[ind]]$mat
      psi <- templates[[ind]]$psi
      n <- runif(1, 1, round(nrow(mat)*ncol(mat))/7)
      mat <- mutate.grid3(mat, n)
    }else if(generate=='neighbors'){
      mat <- templates[[ind]]$mat
      psi <- templates[[ind]]$psi
      asvec <- as.vector(mat)
      maxdist <- sum(sapply(1:length(asvec), function(j){max(abs(5 - asvec[j]), abs(0 - asvec[j]))})^2)^.5
      Dist <- min(maxdist, D)
      mat <- findNeighbors2(mat, Dist, markers)
    }else if(generate=='row-wise'){
      mat <- templates[[ind]]$mat
      psi <- templates[[ind]]$psi
      for(x in 1:nrow(mat)){
        post.curr <- compPost(mat, base.psi, datavec, markers, pars)[1]
        temp <- mat
        sampled <- sample(1:nrow(rowpool), subsample.pool, replace=FALSE)
        post.temps <- rep(NA, length(sampled))
        for(y in 1:length(sampled)){
          temp[x,] <- rowpool[sampled[y],]
          post.temps[y] <- compPost(temp, base.psi, datavec, markers, pars)[1]
        }
        rowpick <- rowpool[sampled[which.max(post.temps)],]
        if(max(post.temps>post.curr)){
          mat[x,] <- rowpick
        }
      }
    }else if(generate=='entry-wise'){
      mat <- templates[[ind]]$mat
      psi <- templates[[ind]]$psi
      asvec <- as.vector(mat)
      subsample.entry <- min(cnmax, subsample.pool)
      for(x in 1:length(asvec)){
        post.curr <- compPost(mat, base.psi, datavec, markers, pars)[1]
        curr <- asvec[x]
        pool <- 0:cnmax[which(!0:cnmax %in% curr)]
        sampled <- sample(pool, subsample.entry, replace=FALSE)
        post.temps <- rep(NA, length(sampled))
        for(y in 1:length(sampled)){
          temp <- asvec
          temp[x] <- pool[sampled[y]]
          asmat <- matrix(temp, nrow=nrow(parset[[ind]]$mat), ncol=ncol(parset[[ind]]$mat))
          post.temps[y] <- compPost(asmat, base.psi, datavec, markers, pars)[1]
        }
        entrypick <- pool[sampled[which.max(post.temps)],]
        if(max(post.temps>post.curr)){
          asvec[x] <- entrypick
        }
      }
      mat <- matrix(asvec, nrow=nrow(parset[[ind]]$mat), ncol=ncol(parset[[ind]]$mat))
    }
    newpost <- compPost(mat, base.psi, datavec, markers, pars)[1]
    keep <- sample(c(TRUE, FALSE), 1, prob=c(p, 1-p))
    if(newpost<curr.post & keep==TRUE){
      mat <- base.mat
    }else{
      curr.post <- newpost
    }
    if(psigen=='findPsi'){
      psi <- try(findPsi(data, mat), silent=TRUE)
      if(class(psi)=='try-error'){
        psis <- rdirichlet(1, alpha=alpha)
      }
    }else if(psigen=='sample'){
      psis <- rdirichlet(subsample.psi, alpha=alpha)
    }else if(psigen=='mutate'){
      psis <- t(sapply(1:subsample.psi, function(j){mutate.psi3(base.psi, deltarange)}))
    }
    psi.res <- t(sapply(1:nrow(psis), function(j){compPost(mat, psis[j,], datavec, markers, pars)}))
    psi.post <- psi.res[,1]
    maxindex <- which.max(psi.post)
    psi <- psis[maxindex,] 
    keep2 <- sample(c(TRUE, FALSE), 1, prob=c(p, 1-p))
    if(psi.post[maxindex]<curr.post & keep2==TRUE){
      psi <- base.psi
    }
    list('mat'=mat, 'psi'=psi)
  })
  out
}

#Actually fitting populations of parameters to distributions:
#...

###Likelihood functions for integrated model:
#CN data
likely1 <- function(model, data, markers, sigma0, log=TRUE){
  cn <- model$cn
  if(class(cn)=='matrix'){
    rsq <- (colSums(model$psi*t(cn)) - data)^2
  }else{
    rsq <- ((model$psi*cn) - data)^2
  }
  sum(dnorm(rsq, 0, sigma0/markers, log=log))
}

#Mutation data
likely2 <- function(model, data, sigma, mu, log=TRUE){
  if(length(model$psi)==1 & sum(model$psi)<1){
    model$psi <- c(1 - model$psi, model$psi)
  }else if(sum(model$psi)<1){
    model$psi <- sort(c(model$psi, 1 - sum(model$psi)), decreasing=TRUE)
  }
  mu.mut <- colSums(model$psi*t(model$cn.mut))*mu
  mu.ref <- colSums(model$psi*t(model$cn.ref))*mu
  rsq <- c(abs(mu.mut - data$varCount), (mu.ref - data$refCount)^2)
  sum(dnorm(rsq, 0, sigma, log=log))
}

#both
likely3 <- function(model, cndata, mutdata, sigma0.seg, sigma.m, mu, log=TRUE){
  cn <- model$cn
  psi <- model$psi
  mut.array <- model$mut.array
  if(class(cn)=='matrix'){
    rsq.cn <- (colSums(model$psi*t(cn)) - cndata)^2
  }else{
    rsq.cn <- ((model$psi*cn) - cndata)^2
  }
  cn.mut <- mut.array[,,1]
  cn.ref <- mut.array[,,2]
  mu.mut <- colSums(model$psi*t(cn.mut))*mu
  mu.ref <- colSums(model$psi*t(cn.ref))*mu
  rsq.m <- c((mu.mut - data$varCount)^2, (mu.ref - data$refCount)^2)
  sum(dnorm(rsq.cn, 0, sigma0.seg/markers, log=log) + 
        dnorm(rsq.m, 0, sigma.m, log=log))
}

#Priors for full integrated model:
prior1 <- function(model, theta, ktheta, stheta, thresh){
  cnPrior(model$cn, model$psi, theta) + psiPrior(model$psi) + kPrior(length(which(model$psi>thresh)), ktheta) +
    sPrior(model$cn, stheta)
}

prior2 <- function(model, theta, ktheta, thresh){
  if(sum(model$psi)<1){
    model$psi <- c(model$psi, 1 - sum(model$psi))
  }
  psiPrior(model$psi) + kPrior(length(which(model$psi>thresh)), ktheta) + 
    cnPrior(model$cn[,,1], model$psi, theta) + cnPrior(model$cn[,,2], model$psi, theta)
}

prior3 <- function(model, theta, ktheta, stheta, thresh){
  cnPrior(model$cn, model$psi, theta) + psiPrior(model$psi) + kPrior(length(which(model$psi>thresh)), ktheta) + 
    sPrior(model$cn, stheta)
}

postComp1 <- function(model, data, markers, pars, log=TRUE){
  sigma0 <- pars$sigma0
  theta <- pars$theta
  ktheta <- pars$ktheta
  stheta <- pars$stheta
  thresh <- pars$thresh
  likely1(model, data, markers, sigma0, log=TRUE) + 
    prior1(model, theta, ktheta, stheta, thresh)
}

postComp2 <- function(model, data, pars, log=TRUE){
  sigma <- pars$sigma
  theta <- pars$theta
  mu <- pars$mu
  likely2(model, data, sigma, mu, log=TRUE) + 
    prior2(model, theta)
}

postComp3 <- function(model, cndata, mutdata, markers, pars, log=TRUE){
  sigma0.seg <- pars$sigma0.seg
  sigma.m <- pars$sigma.m
  mu <- pars$mu
  theta <- pars$theta
  ktheta <- pars$ktheta
  stheta <- pars$stheta
  thresh <- pars$thresh
  likely3(model, cndata, mutdata, sigma0.seg, sigma.m, mu, log=TRUE) + 
    prior3(model, theta, ktheta, stheta, thresh)
}

##Row-wise likelihood for mutations:
mutLike <- function(muts, cndata, mutmodel, cnmodel, psi, pars){
  r1 <- sum(mutmodel$ref*psi)*pars$cov - muts$refCounts
  r2 <- sum(mutmodel$var*psi)*pars$cov - muts$varCounts
  rcn <- sum(mutmodel$ref*psi + mutmodel$var*psi) - cndata
  dnorm(r1, 0, pars$sigma.mut, log=TRUE) + dnorm(r2, 0, pars$sigma.mut, log=TRUE) + 
    dnorm(rcn, 0, pars$sigma.cn, log=TRUE)
}

##Row-wise likelihood for CN data with (possible) mutations
rowLike <- function(cndata, mutdata, cnmodel, psi, pars){
  submodels <- lapply(1:length(cnmodel), function(i){
    ref <- 0:cnmodel[i]
    var <- cnmodel[i] - ref
    cbind(ref, var)
  })
  indextab <- expand.grid(lapply(1:length(submodels),function(j){
    1:nrow(submodels[[j]])}))
  models <- array(NA, dim=c(nrow(indextab),length(cnmodel),2))
  for(i in 1:dim(models)[1]){
    for(j in 1:ncol(indextab)){
      models[i,j,] <- submodels[[j]][indextab[i,j],]
    }
  }
  res <- lapply(1:nrow(mutdata), function(i){
    muts <- mutdata[i,]
    lls <- sapply(1:dim(models)[1], function(j){
      mutmodel <- as.data.frame(models[j,,])
      colnames(mutmodel) <- c('ref', 'var')
      mutLike(muts, cndata, mutmodel, cnmodel, psi, pars)})
    list('ll'=max(lls), 'model'=models[which.max(lls),,])
  })
  picks <- lapply(1:length(res), function(i){res[[i]][[2]]})
  model <- array(NA, dim=c(nrow(mutdata), length(cnmodel), 2))
  for(i in 1:dim(model)[1]){
    model[i,,] <- picks[[i]]
  }
  ll <- sum(sapply(1:length(res), function(i){res[[i]][[1]]}))
  output <- list('model'=model, 'll'=ll)
  output
}

#More generic functions:
rowMax <- function(data) sapply(1:nrow(data),function(i){max(data[i,])})
colMax <- function(data) sapply(1:ncol(data),function(i){max(data[,i])})
rowMin <- function(data) sapply(1:nrow(data),function(i){min(data[i,])})
colMin <- function(data) sapply(1:ncol(data),function(i){min(data[,i])})
runif.uniq <- function(n, a, b){
  pass <- FALSE
  while(pass==FALSE) {
    set <- round(runif(n, a, b))
    pass <- length(unique(set))==n
  }
  set
}
  
###DeepCNV functions:
findGeneMuts <- function(gene, genetab, mutdata){
  index <- which(genetab$name==gene)[1]
  start <- genetab$start[index]
  end <- genetab$end[index]
  chr <- unique(as.character(genetab$chrom[index]))
  subset <- mutdata[which(as.character(mutdata$chr)==chr),]
  df <- subset[which(subset$begin >= start & subset$begin <= end),]
  if(nrow(df)>0){
    df$gene <- rep(gene, nrow(df))
  }
  output <- data.frame('refCount'=df$refCount, 'varCount'=df$varCount, 'status'=df$status,
        'begin'=df$begin, 'chr'=df$chr, 'totalCount'=df$totalCount)
  output
}

###Clustering:
ClusterNum<-function(X, maxh,d="euclidean",method="ward.D",hclust=NULL){
  #X is required
  if(class(hclust)=="NULL"){
    X<-data.frame(scale(X));
    d<-dist(X,d);
    hc<-hclust(d,method);
    N<-min(c(maxh, nrow(X)));
  }else if(class(hclust)=="hclust"){
    hc<-hclust;
    N<-length(hc$order);
  }else{
    stop("The class of hclust is wrong!");
  }  
  R2<-seq(0,1,length=N);
  PSF<-rep(0,N);PST2<-PSF;
  T<-sum(diag(cov(X)*(N-1)));
  Nm<-rep(0,N-1);Wm<-rep(0,N-1);
  pdf("cluster");plot(1,type="n");
  rect0<-rect.hclust(hc,k=2);
  for(t in 2:(N-1)){
    rect<-rect0;
    if(t<N-1){
      rect0<-rect.hclust(hc,k=t+1);
    }
    b<-colMeans(X[rect[[1]],])-colMeans(X);
    B<-length(rect[[1]])*t(b)%*%b;
    for(j in 1:t){
      if(length(rect[[j]])!=length(rect0[[j]])){
        break;}
    }
    Nm[t]<-length(rect[[j]]);
    Wm[t]<-sum(diag(cov(X[rect[[j]],])*(length(rect[[j]])-1)));
    for(i in t:2){
      b<-colMeans(X[rect[[i]],])-colMeans(X);
      B<-B+length(rect[[i]])*t(b)%*%b;
    }
    R2[t]<-B/T;
    PSF[t]<-(B*(N-t))/((T-B)*(t-1));
  }
  dev.off();
  R1<-R2[-1];SPRSQ<-matrix(R1-R2[1:(N-1)]);
  Bkl<-SPRSQ*T;
  PST2<-(Nm-2)*Bkl/(Wm-Bkl);
  PST2[1]<-(N-2)*Bkl[1]/(T-Bkl[1]);
  PST2[c(N-1,N-2)]<-0;
  data.frame(RSQ=R2[1:(N-1)],SPRSQ,PSF=PSF[1:(N-1)],PST2);
}

#Breakpoint detection with variant calls:
seqSegment <- function(seqdata, cutoff=.1){
  sigma <- (sd(seqdata$refCounts) + sd(seqdata$varCounts))/2
  mu <- (mean(seqdata$refCounts) + mean(seqdata$varCounts))/2
  grain <- 100
  segmented <- lapply(unique(seqdata$chr), function(chr){
    data <- seqdata[seqdata$chr==chr,]
    remainder <- nrow(data) - grain*(floor(nrow(data)/grain))
    medians <- t(sapply(1:(floor(nrow(data)/grain)+1),function(i){
      c(median(data$refCounts[(50*(i-1)+1):(50*i)]), median(data$varCounts[(50*(i-1)+1):(50*i)]))
    }))
    colnames(medians) <- c('ref.med', 'var.med')
    diffs.ref <- diff(medians[,1])
    diffs.var <- diff(medians[,1])
    if(remainder>0){
      sigmas <- c(rep(sigma/grain^.5, floor(nrow(data)/grain)-1), sigma/remainder^.5)
    }else{
      sigmas <- c(rep(sigma/grain^.5, floor(nrow(data)/grain)-1))
    }
    p1 <- 2*(1 - pnorm(abs(diffs.ref), 0, sigmas))
    p2 <- 2*(1 - pnorm(abs(diffs.var), 0, sigmas))
    p <- sapply(1:nrow(medians),function(j){min(c(p1[j], p2[j]))})
    segstarts <- c(1, which(p<=cutoff)+1)
    segends <- c(segstarts[-1]-1, nrow(medians))
    starts <- (segstarts - 1)*grain + 1
    ends <- c(starts[-1]-1, (nrow(medians)-1)*grain + remainder)
    markers <- ends - starts + 1
    X <- sapply(1:length(segstarts), function(i){mean(medians[(segstarts[i]:segends[i]),1])/mu})
    Y <- sapply(1:length(segstarts), function(i){mean(medians[(segstarts[i]:segends[i]),2])/mu})
    data.frame('chr'=rep(chr, length(markers)), 'seg'=1:length(markers), 'start'=starts, 
      'end'=ends, 'LRR'=log10((X+Y)/2), 'BAF'=X/(X+Y), 'X'=X, 'Y'=Y, 'markers'=markers)
  })
  segmented <- Reduce(rbind, segmented)
  segmented
}

seqSeg <- function(seqsnps, len, thresh){
  sigma0 <- sd(c(seqsnps$refCounts,seqsnps$varCounts))
  dens <- density(seqsnps$refCounts)
  mu <- dens$x[which.max(dens$y)]
  chrs <- unique(seqsnps$chr)
  dflist <- lapply(1:length(chrs),function(i){
    seqdata.chr <- seqsnps[seqsnps$chr==chrs[i],]
    meandiffs <- t(sapply(1:(nrow(seqdata.chr)-len),function(j){
      mu1.ref <- mean(seqdata.chr$refCounts[((j-1)*len+1):(j*len)])
      mu2.ref <- mean(seqdata.chr$refCounts[(j*len+1):((j+1)*len)])
      mu1.var <- mean(seqdata.chr$varCounts[((j-1)*len+1):(j*len)])
      mu2.var <- mean(seqdata.chr$varCounts[(j*len+1):((j+1)*len)])
      c(abs(mu2.ref - mu1.ref), abs(mu2.var - mu1.var))
    }))
    peaks.ref <- findPeaks(meandiffs[,1], thresh=.25)
    peaks.var <- findPeaks(meandiffs[,2], thresh=.25)
    peaks <- sort(unique(c(peaks.ref, peaks.var)), decreasing=FALSE)
    if(length(peaks)>0){
      start <- c(1, 1+len*sapply(1:length(peaks),function(j){peaks[j]}))
      end <- c(start[-1]-1, nrow(seqdata.chr))
    }else{
      start <- 1
      end <- nrow(seqdata.chr)
    }
    markers <- end - start + 1
    X <- sapply(1:length(start),function(j){median(seqdata.chr$refCounts[start[j]:end[j]])})/mu
    Y <- sapply(1:length(start),function(j){median(seqdata.chr$varCounts[start[j]:end[j]])})/mu
    LRR <- log10((X+Y)/2)
    BAF <- X/(X+Y)
    data.frame('chr'=rep(chrs[i],length(start)), 'start'=start, 'end'=end, 'LRR'=LRR, 
               'BAF'=BAF, 'X'=X, 'Y'=Y, 'markers'=markers)
  })
  df <- Reduce(rbind, dflist)
  df$seg <- 1:nrow(df)
  df
}
