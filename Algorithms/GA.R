source('C:\\Users\\Mark\\OneDrive - The Ohio State University\\Rcode\\functions.r')
source('C:\\Users\\Mark\\OneDrive - The Ohio State University\\Rcode\\loading.r')

bin <- TRUE
if(bin){
  sigma0 <- 1
  data <- datamat
  markers <- rep(datamat$markers, 2)
  theta <- .99
  alpha <- .2
  iters <- 500
  P <- 50
  O <- 5
  M <- 10
  ktheta <- .99
  kmax <- 5
}

#after grid.pop and psi.pop have been defined:
truemat.t <- truemat2
truepsi.t <- c(truepsi, rep(0, kmax - length(truepsi)))
psis.filt <- rbind(psis.filt, c(truepsi, rep(0, kmax-length(truepsi))))
ind.true <- nrow(psis.filt)
dims <- dim(grid.filt)
dims[1] <- dims[1] + 1
temp <- array(data=NA, dim=c(dims))
for(i in 1:(dim(temp)[1]-1)){
  for(j in 1:dim(temp)[2]){
    for(k in 1:dim(temp)[3]){
      temp[i,j,k] <- grid.filt[i,j,k]
    }
  }
}
temp[dims[1],,] <- truemat.t
grid.filt <- temp
resids <- rbind(resids, (colSums(truepsi.t*truemat.t)-datavec)^2)

###
GA <- function(sigma0, theta, alpha, iters, P, O, M, ktheta, data, markers){
  t1 <- Sys.time()
  datavec <- c(data$X, data$Y)
  nmarks <- rep(markers, 2)
  filtered <- filter(psipool, vecpool, datavec)
  kmax <- ncol(psipool)
  psis.filt <- filtered[[1]]
  grid.filt <- filtered[[2]]
  resids <- filtered[[3]]
  rsums <- rowSums(resids)
  ll.init <- sapply(1:nrow(psis.filt), function(i){sum(logLike(resids[i,], nmarks, sigma0))})
  psipriors.init <- sapply(1:nrow(psis.filt), function(i){psiPrior(psis.filt[i,], alpha)})
  cnpriors.init <- sapply(1:nrow(psis.filt), function(i){cnPrior(grid.filt[i,,], theta)})
  ks <- sapply(1:nrow(psis.filt), function(i){length(which(psis.filt[i,]>0))})
  kpriors.init <- sapply(ks, kPrior, ktheta)
  priors.init <- psipriors.init + cnpriors.init + kpriors.init
  #priors.init <- sapply(1:nrow(psis.filt), function(i){cnPrior(grid.filt[i,,], theta) + psiPrior(psis.filt[i,], alpha)})
  posts.init <- ll.init + priors.init
  posts <- posts.init
  ord <- order(posts, decreasing=TRUE)
  grid.pop <- grid.filt
  psi.pop <- psis.filt
  size <- P + 2*O*P
  psi.array <- array(data=NA, dim=c(size,ncol(psi.pop),iters))
  grid.array <- array(data=NA, dim=c(size,length(datavec),kmax,iters))
  post.mat <- matrix(data=NA, nrow=size, ncol=iters)
  ll.mat <- matrix(data=NA, nrow=size, ncol=iters)
  exp.array <- array(data=NA, dim=c(size, length(datavec), iters))
  iter <- 1
  while(iter <= iters){
    probs <- convert(posts)
    picked <- matrix(order(probs, decreasing=TRUE)[1:(2*P)], ncol=2, nrow=P)
    parents1 <- lapply(1:nrow(picked), function(i){grid.pop[picked[i,1],,]})
    parents2 <- lapply(1:nrow(picked), function(i){grid.pop[picked[i,2],,]})
    offspring <- crossover2(parents1, parents2, O)
    mutated <- lapply(1:length(offspring), function(i){mutate.grid(offspring[[i]])})
    mutated.psi <- t(sapply(1:length(mutated), function(i){
      trial <- try(findPsi(datavec, mutated[[i]]), silent=TRUE)
      if(class(trial)=='try-error'){
        output <- psi.pop[nrow(picked) + i,]
      }else{
        output <- trial
      }
      as.vector(output)
    }))
    offspring.psi <- t(sapply(1:length(offspring), function(i){
      trial <- try(findPsi(datavec, offspring[[i]]), silent=TRUE)
      if(class(trial)=='try-error'){
        output <- psi.pop[nrow(picked) + length(mutated) + i,]
      }else{
        output <- trial
      }
      if(length(which(is.na(output)))>0){
        output <- rdirichlet(1, rep(1, kmax))
      }
      as.vector(output)
    }))
    best <- order(posts, decreasing=TRUE)[1:P]
    best.mats <- lapply(1:length(best), function(i){grid.pop[best[i],,]})
    best.psis <- t(sapply(1:length(best), function(i){psipool[best[i],]}))
    newmats <- c(best.mats, offspring, mutated)
    newpsis <- rbind(mutated.psi, offspring.psi)
    newpsis <- rbind(best.psis, newpsis)
    grid.pop <- array(data=NA, dim=c(length(newmats), dim(grid.pop)[2:3]))
    psi.pop <- matrix(data=NA, ncol=kmax, nrow=nrow(newpsis))
    for(i in 1:length(newmats)){
      grid.pop[i,,] <- newmats[[i]]
      psi.pop[i,] <- c(newpsis[i,], rep(0, kmax - length(newpsis[i,])))
    }
    resids <- matrix(data=NA, nrow=nrow(grid.pop), ncol=length(datavec))
    for(i in 1:length(datavec)){
      resids[,i] <- rowSums((grid.pop[,i,] * psi.pop - datavec[i])^2)
    }
    resids <- t(resids)
    lls <- sapply(1:nrow(psi.pop),function(i){sum(logLike(resids[,i], nmarks, sigma0))})
    ks <- sapply(1:nrow(psi.pop), function(i){length(which(psi.pop[i,]>0))})
    psipriors <- sapply(1:nrow(psi.pop), function(i){psiPrior(psi.pop[i,], alpha)})
    cnpriors <- sapply(1:nrow(psi.pop), function(i){cnPrior(grid.pop[i,,], theta)})
    kpriors <- sapply(ks, kPrior, ktheta)
    priors <- psipriors + cnpriors + kpriors
    posts <- lls +  priors
    expmat <- t(sapply(1:size, function(i){rowSums(grid.pop[i,,]*psi.pop[i,])}))
    psi.array[,,iter] <- psi.pop
    grid.array[,,,iter] <- grid.pop
    post.mat[,iter] <- posts
    ll.mat[,iter] <- lls
    exp.array[,,iter] <- expmat
    iter <- iter + 1
    
    #truebest.mat <- findBest.mat(grid.pop, truemat)
    #truebest.psi <- findBest.psi(psi.pop, truepsi)
  }
  list(psi.array, grid.array, post.mat)
  
  maxposts <- sapply(1:ncol(post.mat), function(i){max(post.mat[,i])})
  #plot(maxposts)
  maxlls <- sapply(1:ncol(ll.mat), function(i){max(ll.mat[,i])})
  #plot(maxlls)
  maxindices <- sapply(1:ncol(post.mat), function(i){which.max(post.mat[,i])})
  bestpsi <- psi.array[maxindices[iters],,iters]
  bestmat <- grid.array[maxindices[iters],,,iters]
  bestexp <- exp.array[maxindices[iters],,iters]
  bestpost <- max(maxposts)
  residsums <- sapply(1:iters, function(i){sum(exp.array[maxindices[i],,i])})
  final <- list(bestmat, bestpsi, bestpost)
  arrays <- list(grid.array, psi.array, post.mat, exp.array)
  t2 <- Sys.time()
  tdiff <- t2 - t1
  list('final'=final, 'arrays'=arrays, 'runtime'=tdiff)
}

#res <- GA(sigma0=1, theta=.8, alpha=.5, iters=500, P=50, O=5, M=10, ktheta=.5, data=datamat, markers=markers)

#Note: where I left off: resids, psis.filt, grid.filt are all defined with the last
#vector and matrix corresponding to the truth, which has the 164th highest posterior of
#the original population...while loop has not been started