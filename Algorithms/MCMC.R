library(gtools)
load('C:\\Users\\Mark\\Dropbox\\Lab\\cll\\Simulations\\sims1.rda')
source('C:\\Users\\Mark\\OneDrive - The Ohio State University\\Rcode\\functions.r')
source('C:\\Users\\Mark\\OneDrive - The Ohio State University\\Rcode\\loading.r')
sigma0 <- 1
theta <- .5
ktheta <- .5
alpha <- .5
N <- 300
n <- 50
nmarks <- markers
kmax <- 5
iters <- 100

MCMC <- function(){
  t1 <- Sys.time()
  filtered <- filter(psipool, vecpool, datavec)
  psis.filt <- filtered[[1]]
  psis.pop <- psis.filt
  grid.filt <- filtered[[2]]
  grid.pop <- grid.filt
  resids <- filtered[[3]]
  ll.init <- t(sapply(1:nrow(psis.filt), function(i){logLike(resids[i,], nmarks, sigma0)}))
  psipriors.init <- sapply(1:length(datavec), function(j){sapply(1:nrow(psis.filt), function(i){psiPrior(psis.filt[i,], alpha)})})
  cnpriors.init <- sapply(1:length(datavec), function(j){sapply(1:nrow(psis.filt), function(i){cnPrior(grid.filt[i,j,], theta)})})
  ks <- sapply(1:nrow(psis.filt), function(i){length(which(psis.filt[i,]>0))})
  kpriors.init <- sapply(1:length(datavec), function(j){sapply(ks, kPrior, ktheta)})
  priors.init <- psipriors.init + cnpriors.init + kpriors.init
  #priors.init <- sapply(1:nrow(psis.filt), function(i){cnPrior(grid.filt[i,,], theta) + psiPrior(psis.filt[i,], alpha)})
  posts.init <- ll.init + priors.init
  posts.orig <- posts.init
  psi.array <- array(data=NA, dim=c(length(datavec), kmax, N*iters))
  grid.array <- array(data=NA, dim=c(length(datavec), kmax, N*iters))
  post.mat <- matrix(data=NA, ncol=length(datavec), nrow=N*iters)
  vecs <- array(data=NA, dim=c(N,length(datavec),kmax))
  psis <- array(data=NA, dim=c(N,length(datavec),kmax))
  posts <- matrix(data=NA, nrow=N, ncol=(length(datavec)))
  iter <- 1
  while(iter <= iters) {
    curr.max <- N*iter
    curr.min <- N*(iter-1)+1
    for(j in 1:length(datavec)){
      if(iter==1){
        ord <- order(posts.orig[,j], decreasing=TRUE)
        indices <- ord[1:N]
        psipop <- psis.pop[c(indices),]
        gridpop <- grid.pop[c(indices),j,] 
        post.mat[curr.min:curr.max,j] <- posts.orig[indices,j]
      }else{
        psipop <- psis[,j,]
        gridpop <- vecs[,j,]
        post.mat[curr.min:curr.max,j] <- posts[,j]
      }
      psi.array[j,,curr.min:curr.max] <- t(psipop)
      grid.array[j,,curr.min:curr.max] <- t(gridpop)
      newvecs <- resample.grid(t(grid.array[j,,1:(curr.max)]), post.mat[1:(curr.max),j], N)
      newpsis <- resample.psi(t(psi.array[j,,1:(curr.max)]), post.mat[1:(curr.max),j], N)
      for(i in 1:nrow(newvecs)){
        newvecs[i,c(which(newpsis[i,]==0))] <- rep(0, length(which(newpsis[i,]==0)))
      }
      newresids <- rowSums(newvecs*newpsis) - datavec[j]
      newlikes <- logLike(newresids, nmarks[j], sigma0)
      newpriors <- sapply(1:nrow(newpsis), function(k){psiPrior(newpsis[k,], alpha) + cnPrior(newvecs, theta) +
                                                         kPrior(length(which(newpsis[k,]>0)), ktheta)})
      newposts <- newlikes + newpriors
      vecs[,j,] <- newvecs
      psis[,j,] <- newpsis
      posts[,j] <- newposts
    }
    iter <- iter + 1
  }
  arrays <- list(psi.array, grid.array, post.mat)
  postsums <- colSums(post.mat)
  #smoothed <- (postsums)
  max.index <- which.max(postsums)
  maxpost <- max(postsums)
  bestpsi <- psi.array[max.index,,]
  bestmat <- grid.array[max.index,,,]
  final <- list(bestmat, bestpsi, maxpost)
  t2 <- Sys.time()
  tdiff <- t2 - t1
  list('final'=final, 'arrays'=arrays, 'runtime'=tdiff)
}
