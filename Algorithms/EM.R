source('C:\\Users\\Mark\\OneDrive - The Ohio State University\\Rcode\\functions.r')
source('C:\\Users\\Mark\\OneDrive - The Ohio State University\\Rcode\\loading.r')
#N <- 300
#X <- 20
#a <- 1
#b <- 100
#n <- 5
#ktheta <- .5
#ntheta <- .01
#theta <- .2
#epsilon <- 10
#max.iters <- 150
#kmax <- 5
#sigma0 <- 1
#alpha <- .5
#data <- datamat

###
EM <- function(N, X, a, b, n, theta, ntheta, ktheta, epsilon, max.iters, kmax, sigma0, alpha, data, markers){
  t1 <- Sys.time()
  datavec <- c(data$X, data$Y)
  nmarks <- rep(markers, 2)
  filtered <- filter(psipool, vecpool, datavec)
  kmax <- ncol(psipool)
  psis.filt <- filtered[[1]]
  grid.filt <- filtered[[2]]
  newgrid <- grid.filt
  newpsis <- psis.filt
  resids <- filtered[[3]]
  ll.init <- sapply(1:nrow(psis.filt), function(i){sum(logLike(resids[i,], nmarks, sigma0))})
  psipriors.init <- sapply(1:nrow(psis.filt), function(i){psiPrior(psis.filt[i,], alpha)})
  cnpriors.init <- sapply(1:nrow(psis.filt), function(i){cnPrior(grid.filt[i,,], theta)})
  ks <- sapply(1:nrow(psis.filt), function(i){length(which(psis.filt[i,]>0))})
  kpriors.init <- sapply(ks, kPrior, ktheta)
  priors.init <- psipriors.init + cnpriors.init + kpriors.init
  #priors.init <- sapply(1:nrow(psis.filt), function(i){cnPrior(grid.filt[i,,], theta) + psiPrior(psis.filt[i,], alpha)})
  posts.init <- ll.init + priors.init
  likes <- ll.init
  priors <- priors.init
  posts.curr <- likes + priors
  bestpost <- max(posts.curr)
  diff <- epsilon + 1
  iter <- 1
  psi.mat <- matrix(data=NA, nrow=max.iters, ncol=kmax)
  mat.array <- array(data=NA, dim=c(max.iters,kmax,length(datavec)))
  postvec <- rep(NA, max.iters)
  diffs <- rep(NA, max.iters)
  expmat <- matrix(data=NA, nrow=max.iters, ncol=length(datavec))
  while(iter <= max.iters) {
    last.bestpost <- bestpost
    bestn <- order(posts.curr, decreasing=TRUE)[1:n]
    #probs <- convert(posts.curr)
    #bases <- sample(1:length(probs), N, replace=TRUE, prob=probs)
    if(iter==1){
      bases <- sample(1:length(posts.curr), X*N, replace=FALSE)
    }else{
      bases <- rep(1:N, X)
    } 
    deltas <- rbeta(N*X, a, b)
    nprobs <- dexp(1:length(grid.filt[1,,]), ntheta)
    nset <- sample(1:length(grid.filt[1,,]), length(bases), replace=TRUE, prob=nprobs)
    gridlist <- lapply(1:length(bases), function(j){mutate.grid3(grid.filt[bases[j],,], nset[j])})
    tempgrid <- array(data=NA, dim=c(N*X+n,kmax,length(datavec)))
    for(i in 1:length(gridlist)){
      tempgrid[i,,] <- gridlist[[i]]
    }
    for(i in 1:n){
      tempgrid[N*X+i,,] <- newgrid[bestn[i],,]
    }
    newgrid <- tempgrid
    temppsis <- t(sapply(1:length(bases), function(j){mutate.psi(newpsis[bases[j],], deltas[j])}))
    newpsis <- rbind(temppsis, newpsis[c(bestn),])
    resids <- t(sapply(1:nrow(newpsis), function(i){colSums(newgrid[i,,]*newpsis[i,])}))
    newlikes <- sapply(1:nrow(newpsis), function(i){sum(logLike(resids[i,], nmarks, sigma0))})
    psipriors <- sapply(1:nrow(newpsis), function(i){psiPrior(newpsis[i,], alpha)})
    cnpriors <- sapply(1:nrow(newpsis), function(i){cnPrior(newgrid[i,,], theta)})
    ks <- sapply(1:nrow(newpsis), function(j){length(which(newpsis[j,]>0))})
    kpriors <- sapply(ks, kPrior, ktheta=ktheta)
    newpriors <- psipriors + cnpriors
    posts.curr <- newlikes + newpriors
    index.max <- which.max(posts.curr)
    bestpost <- max(posts.curr)
    bestpsi <- newpsis[index.max,]
    bestmat <- newgrid[index.max,,]
    exp <- colSums(bestmat*bestpsi)
    diff <- bestpost - last.bestpost
    psi.mat[iter,] <- bestpsi
    mat.array[iter,,] <- bestmat
    diffs[iter] <- diff
    postvec[iter] <- bestpost
    expmat[iter,] <- exp
    iter <- iter + 1
  }
  final <- list(bestmat, bestpsi, bestpost)
  arrays <- list(mat.array, psi.mat, postvec, expmat)
  t2 <- Sys.time()
  tdiff <- t2 - t1
  list('final'=final, 'arrays'=arrays, 'runtime'=tdiff)
}

