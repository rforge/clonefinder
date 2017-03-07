source('C:\\Users\\Mark\\OneDrive - The Ohio State University\\Rcode\\functions.r')
source('C:\\Users\\Mark\\OneDrive - The Ohio State University\\Rcode\\loading.r')
ktheta <- .5
ntheta <- .01
theta <- .2
epsilon <- 10
kmax <- 5
sigma0 <- 1
alpha <- .5
max.iters <- 100
starting.delta <- .1
n <- 10
N <- 20

###
Alg5 <- function(theta, ktheta, ntheta, epsilon, kmax, sigma0, alpha, max.iters, starting.delta, n, N, data, markers){
  t1 <- Sys.time()
  datavec <- c(data$X, data$Y)
  nmarks <- rep(markers, 2)
  psipicks <- matrix(data=NA, nrow=n, ncol=kmax)
  postpicks <- rep(NA, n)
  matpicks <- array(data=NA, dim=c(length(datavec),kmax,n))
  filtered <- filter(psipool, vecpool, datavec)
  kmax <- ncol(psipool)
  psis.filt <- filtered[[1]]
  grid.filt <- filtered[[2]]
  resids <- filtered[[3]]
  residsums <- rowSums(resids)
  nbest <- order(residsums)[1:n]
  starting.mats <- grid.filt[c(nbest),,]
  starting.psis <- psis.filt[c(nbest),]
  starting.priors <- sapply(1:n,function(i){cnPrior(starting.mats[i,,], theta) + 
    kPrior(length(which(starting.psis[i,]>0)), ktheta) + psiPrior(starting.psis[i,], alpha)})
  starting.lls <- sapply(1:n,function(i){logLike(resids[nbest[i],], markers, sigma0)})
  starting.posts <- starting.priors + starting.lls
  eligible.start <- expand.grid(c(1:kmax), c(1:length(datavec)))
  eligible.start <- rbind(eligible.start, eligible.start)
  eligible.start$dir <- c(rep(1, nrow(eligible.start)/2), rep(-1, nrow(eligible.start)/2))
  colnames(eligible.start)[1:2] <- c('clone', 'segment')
  psi.eligible <- expand.grid(c(1:kmax), c(-1,1))
  colnames(psi.eligible) <- c('clone', 'dir')
  for(i in 1:n){
    mat.start <- starting.mats[i,,]
    psi.start <- starting.psis[i,]
    delta <- starting.delta
    jmats <- array(data=NA, dim=c(length(datavec),kmax,max.iters))
    jpsis <- matrix(data=NA, nrow=max.iters, ncol=kmax)
    jposts <- rep(NA, max.iters)
    oldbestpost <- starting.posts[i]
    j <- 1
    while(j <= max.iters) {
      kmats <- array(data=NA, dim=c(kmax, length(datavec), nrow(psi.eligible)))
      kposts <- rep(NA, nrow(psi.eligible))
      kpsis <- matrix(data=NA, nrow=nrow(psi.eligible), ncol=kmax)
      for(k in 1:nrow(psi.eligible)){
        eligible <- eligible.start
        psi <- mutate.psi2(psi.start, index=eligible$clone[i], delta=delta, dir=eligible$dir[i])
        resids.start <- colSums((psi*mat.start)^2)
        cloneprior <- psiPrior(psi, alpha) + kPrior(length(which(psi>0)), ktheta)
        prior.start <- cnPrior(mat.start, theta) + cloneprior
        ll.start <- sum(logLike(resids.start, markers, sigma0))
        post.start <- ll.start + prior.start
        lmats <- array(data=NA, dim=c(kmax,length(datavec),N))
        lposts <- rep(NA, N)
        l <- 1
        while(l <= N) {
          mat <- mat.start
          pick <- sample(1:nrow(eligible), 1)
          seg <- eligible$segment[pick]
          sign <- eligible$dir[pick]
          clone <- eligible$clone[pick]
          mat[clone,seg] <- mat[clone,seg] + sign
          resids <- colSums((psi*mat)^2)
          prior <- cnPrior(mat, theta) + psiPrior(psi, alpha) + cloneprior
          ll <- sum(logLike(resids, markers, sigma0))
          post <- ll + prior
          if(post < post.start){
            slated <- which(eligible$segment==seg & eligible$dir==sign & 
                            eligible$clone >= clone)
            retain <- which(!1:nrow(eligible) %in% slated)
            eligible <- eligible[c(retain),]
          }
          lmats[,,l] <- mat
          lposts[l] <- post
          l <- l + 1
        }
        bestpost <- max(lposts)
        index.best <- which.max(lposts)
        bestmat <- lmats[,,index.best]
        kmats[,,k] <- bestmat
        kposts[k] <- bestpost
        kpsis[k,] <- psi
      }
      newbestpost <- max(kposts)
      index.best2 <- which.max(kposts)
      best.psi <- kpsis[index.best2,]
      best.mat <- kmats[,,index.best2]
      if(newbestpost > oldbestpost){
        psi.start <- best.psi
        mat.start <- best.mat
        oldbestpost <- newbestpost
      }else{
        delta <- delta/2
      }
      jmats[,,j] <- best.mat
      jpsis[j,] <- best.psi
      jposts[j] <- newbestpost
      j <- j + 1 
    }
    index.best3 <- which.max(jposts)
    matpicks[,,i] <- jmats[,,index.best3]
    psipicks[i,] <- jpsis[index.best3,]
    postpicks[i] <- jposts[index.best3]
  }
  index <- which.max(postpicks)
  best.post <- postpicks[index]
  best.mat <- matpicks[,,index]
  best.psi <- psipicks[index,]
  #null <- scramble(best.mat, best.psi)
  #nullposts <- sapply(1:length(null), function(){compute.post()})
  #which(nullpost > best.[post])
  final <- list(best.mat, best.psi, best.post)
  arrays <- list(matpicks, psipicks, postpicks)
  t2 <- Sys.time()
  tdiff <- t2 - t1
  list('final'=final, 'arrays'=arrays, 'runtime'=tdiff)
}

res <- Alg5(theta=.2, ktheta=.5, ntheta=.01, epsilon=10, kmax=5, sigma0=1, alpha=.5, max.iters=100, 
            starting.delta=.1, n=10, N=20, data=datamat, markers=markers)
