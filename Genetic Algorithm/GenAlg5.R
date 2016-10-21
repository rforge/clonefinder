#'simulator' is needed to generate the simulations
source('C:\\Users\\Mark\\Dropbox\\R Files\\cll\\simulator.r')
library(gtools)
#Functions:

#Find psi given data and an AB matrix
findPsi <- function(data, AB){
  AB <- as.matrix(AB)
  M <- rbind(AB[,1:(ncol(AB)/2)], AB[,(ncol(AB)/2+1):ncol(AB)])
  dat <- as.matrix(c(data[,1], data[,2]))
  psi <- as.vector(solve(t(M) %*% M) %*% t(M) %*% dat)
  psi <- psi/sum(psi)
  psi
}

#In case we need to generate 'segments' or matrices:
generateSegs <- function(data, i, k, theta, priordf=NULL){
  if(is.null(priordf)){
    probs <- dgeom(abs(0:5 - 1), theta)
    Apicks <- sample(0:5, i*k, prob=probs, replace=TRUE)
    Bpicks <- sample(0:5, i*k, prob=probs, replace=TRUE)
    A <- matrix(Apicks, nrow=i, ncol=k)
    B <- matrix(Bpicks, nrow=i, ncol=k)
    AB <- as.data.frame(cbind(A, B))
    output <- AB
  }else{
    cnmat <- bank[,c((k+3):ncol(priordf))]
    lposts <- bank[,1]
    smat <- sapply(1:ncol(cnmat), function(j){sapply(0:5, function(n){sum(lposts[which(
      cnmat[,j]==n)])})})
    pmat <- exp(smat)
    ok <- which(sapply(1:ncol(pmat), function(j){length(which(pmat==0))==0}))
    for(j in 1:ncol(smat)){
      if(j %in% ok){
        pmat[,j] <- pmat[,j]/sum(pmat[,])
      }else{
        pmat[,j] <- 1/(-smat[,j])/sum(1/(-smat[,j]))
      }
    }
    for(j in 1:ncol(pmat)){
      for(i in 1:nrow(pmat)){
        if(is.na(pmat[i,j])){
          pmat[i,j] <- 1
        }
      }
      pmat[,j] <- pmat[,j]/sum(pmat[,j])
    }
    output <- sapply(1:ncol(smat),function(j){sample(0:5, i, prob=pmat[,j], replace=TRUE)})
    output <- as.data.frame(output)
  }
  output
}

#Log posterior function
findPost <- function(datavec, marks, psi, paramvec, sigma0, theta){
  k <- length(paramvec)/2
  sigma <- sigma0/(markers^.5)
  A <- unlist(paramvec[1:k])
  B <- unlist(paramvec[(k+1):(2*k)])
  exp <- c(sum(psi*A), sum(psi*B))
  ll <- dnorm(unlist(exp - datavec), 0, sigma, log=TRUE)
  lprior.psi <- log(ddirichlet(psi, rep(1/2, length(psi))))
  lprior.ab <- log(dgeom(c(abs(1-A), abs(1-B)), theta))
  lpost <- sum(ll, lprior.psi, lprior.ab)
  lpost
}

#crossover operator.
crossover <- function(mat1, mat2, mut.rate, theta, m){
  mat1 <- as.matrix(mat1)
  mat2 <- as.matrix(mat2)
  from1 <- 1:(sample(1:(nrow(mat1)-1), 1))
  from2 <- (from1[length(from1)]+1):nrow(mat1)
  probvec <- dgeom(0:m, theta)
  output1 <- matrix(rep(NA, nrow(mat1)*ncol(mat1)), nrow=nrow(mat1), ncol=ncol(mat1))
  output2 <- matrix(rep(NA, nrow(mat1)*ncol(mat1)), nrow=nrow(mat1), ncol=ncol(mat1))
  for(j in from1){
    output1[j,] <- mat1[j,]
    output2[j,] <- mat2[j,]
    toMutate1 <- sample(c(TRUE, FALSE), 1, prob=c(mut.rate, 1-mut.rate))
    toMutate2 <- sample(c(TRUE, FALSE), 1, prob=c(mut.rate, 1-mut.rate))
    if(toMutate1){
      output1[j,] <- sample(0:m, ncol(mat1), prob=probvec, replace=TRUE)
    }
    if(toMutate2){
      output2[j,] <- sample(0:m, ncol(mat1), prob=probvec, replace=TRUE)
    }
  }
  for(j in from2){
    output1[j,] <- mat2[j,]
    output2[j,] <- mat1[j,]
    toMutate1 <- sample(c(TRUE, FALSE), 1, prob=c(mut.rate, 1-mut.rate))
    toMutate2 <- sample(c(TRUE, FALSE), 1, prob=c(mut.rate, 1-mut.rate))
    if(toMutate1){
      output1[j,] <- sample(0:m, ncol(mat1), prob=probvec, replace=TRUE)
    }
    if(toMutate2){
      output2[j,] <- sample(0:m, ncol(mat1), prob=probvec, replace=TRUE)
    }
  }
  list(output1, output2)
}

#Now we generate a simulated tumor, with concomitant 'data' and 'truth'
truepsi <- c(.75, .25)
test1 <- simulate(psilist=list(truepsi), theta=1, sigma=.1, min=.005, snprate=.6)
data <- test1[[3]]
fulldata <- test1[[8]]
#fulldata$pos.gen <- sapply(1:nrow(fulldata), function(j){sum(sapply(1:((fulldata[j,]$chr.n-1)), 
#                    function(i){max(fulldata[fulldata$chr.n==i,]$pos)})) + fulldata$pos[j]})
plot(fulldata[fulldata$chr.n==2,]$pos, fulldata[fulldata$chr.n==2,]$lrr, pch='.', 
     xlab='pos', ylab='lrr')
plot(fulldata[fulldata$chr.n==2,]$pos, fulldata[fulldata$chr.n==2,]$baf, pch='.', 
     xlab='pos', ylab='baf')
plot(fulldata[fulldata$chr.n==2,]$pos, fulldata[fulldata$chr.n==2,]$X, pch='.', 
     xlab='pos', ylab='X')
plot(fulldata[fulldata$chr.n==2,]$pos, fulldata[fulldata$chr.n==2,]$Y, pch='.', 
     xlab='pos', ylab='Y')

truemat <- test1[[5]]
nseg <- nrow(data)
markers <- data$markers
colorset <- c('cadetblue', 'cadetblue1', 'chartreuse1', 'chartreuse4', 'coral', 'cornflowerblue', 
  'cyan', 'darkblue', 'darkgoldenrod','darkgoldenrod1', 'brown1', 'darkgreen', 'burlywood1', 
  'dimgray', 'darkorchid2', 'black', 'yellow', 'firebrick2', 'darkolivegreen1', 'khaki')

iter <- 1
iters <- 100
#k is number of subclones
k <- length(truepsi)
#N is the number of matrices in the 'population'
N <- 500
N.init <- 10000
#q is the number of offspring generated per iteration; it must be an even number, so two can be generated per crossover
q <- 450
theta <- .1
sigma0 <- 10
M <- 5
matlist <-lapply(1:N.init, function(iter){generateSegs(data, nseg, k, theta, priordf=NULL)})
#psis <- t(sapply(1:length(mats), function(j){
#  t <- try(findPsi(data, as.matrix(mats[[j]])), silent=T)
#  if(class(t)=='try-error'){
#    psi <- rdirichlet(1, rep(1/k, k))
#  }else{
#    psi <- abs(t)/sum(abs(t))
#  }
#}))
#posts <- t(sapply(1:nrow(psis), function(j){sapply(1:nrow(data), function(i){findPost(datavec=
#    unlist(data[i,1:2]),marks=data$markers[i], psi=psis[j,], paramvec=mats[[j]][i,], sigma0, theta)})}))
#cnarray stores copy number matrices for each matrix in the population, for each iteration
cnarray <- array(data=NA, dim=c(nrow(data),2*k,N,iters))
mut.rate <- .1
accs <- matrix(rep(NA, iters*N), nrow=iters, ncol=N)
seg.acc.array <- array(data=NA, dim=c(nrow(data),N,iters))
psiarray <- array(data=NA, dim=c(k,N,iters))
exparray <- array(data=NA, dim=c(nrow(data),2,N,iters))
lparray <- array(data=NA, dim=c(nrow(data),N,iters))
lpsums <- matrix(rep(NA, N*iters), nrow=iters, ncol=N)

while(iter <= iters) {
  if(iter==1){
    time1 <- Sys.time()
  }
  logpostmat <- t(sapply(1:length(matlist), function(j){
    t2 <- try(findPsi(data[,1:2], matlist[[j]]), silent=T)
    if(class(t2)=='try-error'){
      imputedPsi <- as.vector(rdirichlet(1, rep(1/2, k)))
    }else{
      imputedPsi <- t2
    }
    temp <- sapply(1:nrow(data), function(i){
      post <- findPost(datavec=unlist(data[i,1:2]), marks=markers[i], psi=imputedPsi, 
               paramvec=matlist[[j]][i,], sigma0, theta)
      post
      }
    )
    c(temp, imputedPsi)
  }))
  logposts <- logpostmat[,1:nrow(data)]
  psis <- logpostmat[,(nrow(data)+1):ncol(logpostmat)]
  pre.curr.lpsums <- rowSums(logposts)
  b <- N - q
  bestb <- order(pre.curr.lpsums, decreasing=TRUE)[1:b]
  curr.lpsums <- pre.curr.lpsums[bestb]
  probs <- pre.curr.lpsums
  maxim <- max(probs)
  probs <- probs - maxim
  probs <- exp(probs)
  minim <- min(probs)
  probs[which(probs==0)] <- minim*(1/sum(probs) + length(which(probs==0))*minim)
  probs <- probs/sum(probs)
  parents <- t(sapply(1:length(matlist), function(j){sample(1:length(probs), 2, replace=FALSE, prob=probs)}))
  offspring <- lapply(1:(q/2), function(j){crossover(matlist[[parents[j,1]]],
                                                 matlist[[parents[j,2]]], mut.rate=mut.rate, theta=.5, m=M)})
  offspring <- unlist(offspring, recursive=FALSE)
  newmatset <- c(matlist[bestb], offspring)
  if(iter > 1){
    lpsums[iter,] <- pre.curr.lpsums
    for(j in 1:N){
      cnarray[,,j,iter] <- as.matrix(matlist[[j]])
      psiarray[,j,iter] <- psis[j,]
      accs[iter,j] <- sum((matlist[[j]] - truemat)^2)
      seg.acc.array[,j,iter] <- rowSums((matlist[[j]] - truemat)^2)
      exparray[,,j,iter] <- cbind(colSums(t(psis[j,]) %*% t(as.matrix(matlist[[j]][,1:k]))), 
                            colSums(t(psis[j,]) %*% t(as.matrix(matlist[[j]][,(k+1):(2*k)]))))
      lparray[,j,iter] <- logposts[j,]
    }
  }
  matlist <- newmatset
  iter <- iter + 1
  if(iter==iters){
    time2 <- Sys.time()
    runtime <- time2 - time1
  }
}

#res <- list()
#save(res, file='C:\\Users\\Mark\\Dropbox\\R Files\\cll\\res.rda')
#res <- sapply(1:nrow(pardf), foo)

#Looking at results:
best <- which.max(lpsums[iters,])
bestmat <- cnarray[,,best,iters]
bestpsi <- psiarray[,best,iters]
trueposts <- sapply(1:nrow(data), function(j){findPost(data[j,1:2], data$markers[j], 
                               findPsi(data, truemat), truemat[j,], sigma0, theta)})
truepost.sum <- sum(trueposts)
#accs[iters,]
#lpsums[iters,]

lpmaxes <- sapply(1:iters, function(j){max(lpsums[j,])})
lpmeds <- sapply(1:iters, function(j){median(lpsums[j,])})
lpmins <- sapply(1:iters, function(j){min(lpsums[j,])})
lpdiffs <- sapply(1:iters, function(j){max(lpsums[j,] - median(lpsums[j,]))})
map.indices <- unlist(sapply(1:iters, function(j){which.max(lpsums[j,])}))
resids <- sapply(1:iters, function(j){sum((t(exparray[,,map.indices[j],j]) - data[,1:2])^2)})
plot(lpmaxes, main='Max Log Posteriors', xlab='Iteration', ylab='Max Log Posterior')
plot(lpmeds, main='Median Log Posteriors', xlab='Iteration', ylab='Median Log Posterior')
#plot(lpmins, main='Min Log Posteriors', xlab='Iteration', ylab='Min Log Posterior')
plot(lpdiffs, main='Max - Med. Log Posteriors', xlab='Iteration', ylab='Max - Med. Log Posterior')
plot(resids, main='Sum of Squared Residuals' , xlab='Iteration', ylab='Sum of R^2')
#plotfun1 <- function(i){
#  best.acc.row <- sapply(1:iters, function(j){seg.acc.array[i,which.max(lpsums[j,]),j]})
#  plot(best.acc.row, main=paste('Row ', i, ' Accuracy', sep=''), xlab='Iteration', 
#       ylab='Row Accuracy', ylim=c(0, 1))
#}
#sapply(1:nrow(data), plotfun1)

#Pvals:
trueorder <- match(bestpsi, sort(bestpsi, decreasing=TRUE))
reordered <- bestmat[,c(trueorder, trueorder + k)]

#m is the max copy number; n is the number of iterations
computePval <- function(truemat, inferredmat, theta, n, m){
  k <- ncol(truemat)/2
  i <- nrow(truemat)
  probs <- dgeom(abs(0:m - 1), theta)
  iterate <- function(j){
    Apicks <- sample(0:m, i*k, prob=probs, replace=TRUE)
    Bpicks <- sample(0:m, i*k, prob=probs, replace=TRUE)
    A <- matrix(Apicks, nrow=i, ncol=k)
    B <- matrix(Bpicks, nrow=i, ncol=k)
    AB <- as.matrix(cbind(A, B))
    sum((AB - truemat)^2)
  }
  sumr2s.null <- sapply(1:n, iterate)
  sumr2 <- sum((inferredmat - truemat)^2)
  1-length(which(sumr2s.null>=sumr2))/n
}

computePval.seg <- function(truevec, inferredvec, theta, n, m){
  k <- length(truevec)/2
  probs <- dgeom(abs(0:m - 1), theta)
  nullsegs <- t(sapply(1:n, function(q){sample(0:m, 2*k, replace=TRUE, prob=probs)}))
  sumr2s.null <- sapply(1:nrow(nullsegs), function(j){sum((nullsegs[j,] - truevec)^2)})
  sumr2 <- sum((inferredvec - truevec)^2)
  1-length(which(sumr2s.null>=sumr2))/n
}

pval.mat <- computePval(truemat, bestmat, theta, 1000, 5)
pvals.segs <- sapply(1:nrow(truemat), function(x){computePval.seg(
  truevec=truemat[x,], inferredvec=bestmat[x,], theta=theta, n=1000, m=5)})
