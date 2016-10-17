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

#mutation operator
mutate <- function(mat, theta=.2, max=5){
  row <- sample(1:nrow(mat), 1)
  probs <- dgeom(abs(1 - 0:max), theta)
  newrow <- sample(0:max, ncol(mat), prob=probs, replace=TRUE)
  newmat <- mat
  newmat[row,] <- newrow
  newmat
}


#mutation operator in which mutations are not purely random
mutate2 <- function(mat){
  
}

#crossover operator; currently only crosses ver one row at a time
crossover <- function(mat1, mat2){
  from1 <- sample(1:nrow(mat1), runif(1, 1, (nrow(mat1)-1)))
  from2 <- which(!1:nrow(mat1) %in% from1)
  mat <- matrix(rep(NA, nrow(mat1)*ncol(mat1)), nrow=nrow(mat1), ncol=ncol(mat1))
  for(j in from1){
    mat[j,] <- mat1[j,]
  }
  for(j in from2){
    mat[j,] <- mat2[j,]
  }
  mat
}

#Now we generate a simulated tumor, with concomitant 'data' and 'truth'
truepsi <- c(.7, .3)
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
iters <- 200
#k is number of subclones
k <- 2
#N is the number of matrices in the 'population'
N <- 100
theta <- .5
mats <-lapply(1:N, function(iter){generateSegs(data, nseg, 2, theta, priordf=NULL)})
#cnarray stores copy number matrices for each matrix in the population, for each iteration
cnarray <- array(data=NA, dim=c(nrow(data),2*k,N,iters))
for(i in 1:length(mats)){
  cnarray[,,i,1] <- as.matrix(mats[[i]])
}
sigma0 <- 1
q1 <- 10
q2 <- 5
q3 <- 5
Q <- q1*q2*(q3+1)+N
accs <- matrix(rep(NA, iters*N), nrow=iters, ncol=N)
seg.acc.array <- array(data=NA, dim=c(nrow(data),N,iters))
psiarray <- array(data=NA, dim=c(k,N,iters))
exparray <- array(data=NA, dim=c(2,nrow(data),N,iters))
lparray <- array(data=NA, dim=c(nrow(data),N,iters))
progeny <- array(data=NA, dim=c(nrow(data),2*k,Q,iters))
postbank <- array(data=NA, dim=c(nrow(data),Q,(iters)))
lpsums <- matrix(rep(NA, N*iters), nrow=iters, ncol=N)

while(iter <= iters) {
  if(iter==1){
    time1 <- Sys.time()
  }
  for(j in 1:N){
    t <- try(findPsi(data, cnarray[,,j,iter]), silent=T)
    if(class(t)=='try-error'){
      while(class(t)=='try-error'){
        cnarray[,,j,iter] <- generateSegs(data, nseg, 2, theta)
        t <- try(findPsi(data, cnarray[,,j,iter]), silent=T)
      }
    }
    if(length(which(t<0))>0 | length(which(t>1))>0){
      t <- abs(t)/sum(abs(t))
    }
    psiarray[,j,iter] <- t
  }
  for(j in 1:N){
    exparray[,,j,iter] <- cbind(cnarray[,1:k,j,iter] %*% psiarray[,j,iter], 
                              cnarray[,(k+1):(2*k),j,iter] %*% psiarray[,j,iter])
  }
  if(iter==1){
    lparray[,,iter] <- sapply(1:N,function(i){sapply(1:nrow(data), function(j){findPost(datavec=unlist(data[
      j,1:2]), marks=markers[j], psi=psiarray[,i,iter], paramvec=cnarray[,,i,iter], sigma0, theta)})})
    curr.lpsums <- colSums(lparray[,,iter])
    lpsums[iter,] <- curr.lpsums
  }
  accs[iter,] <- sapply(1:N,function(j){length(which(as.vector(
    cnarray[,,j,iter]==truemat)))/(ncol(cnarray[,,j,iter])*nrow(cnarray[,,j,iter]))})
  seg.acc.array[,,iter] <-  sapply(1:N,function(j){sapply(1:nrow(data), function(j){
    length(which(cnarray[,,,iter]==truemat[j,]))/(2*k)})})
  current.lps <- sapply(1:N,function(j){sapply(1:nrow(data),function(i){findPost(datavec=data[i,1:2], 
            marks=markers[i], psi=psiarray[,j,iter], paramvec=cnarray[i,,j,iter], sigma0, theta)})})
  mat.lps <- colSums(current.lps)
  lps.new <- -mat.lps/max(mat.lps)
  probs <- exp(lps.new)
  probs <- probs/sum(probs)
  picked <- sample(1:N, q1, prob=probs, replace=FALSE)
  matlist <- list(NULL)
  for(j in 1:q1){
    subset <- which(1:N != picked[j])
    mates <- sample(subset, q2, prob=probs[subset])
    for(i in 1:q2){
      offspring <- crossover(cnarray[,,picked[j],iter], cnarray[,,mates[i],iter])
      for(l in 1:(q3+1)){
        if(l<=q3){
          newmat <- mutate(offspring)
        }else{
          newmat <- offspring
        }
      matlist <- c(matlist, list(newmat))
      }
    }
  }
  matlist <- matlist[2:length(matlist)]
  matlist <- c(matlist, lapply(1:N, function(j){cnarray[,,j,iter]}))
  logposts <- sapply(1:length(matlist), function(j){
    sapply(1:nrow(data), function(i){
      t2 <- try(findPsi(data[,1:2], matlist[[j]]), silent=T)
      if(class(t2)=='try-error'){
        imputedPsi <- rdirichlet(1, rep(1/2, k))
      }else{
        imputedPsi <- t2
      }
      findPost(datavec=unlist(data[i,1:2]), marks=markers[i], psi=imputedPsi, 
               paramvec=matlist[[j]][i,], sigma0, theta)
      })
  })
  pre.curr.lpsums <- colSums(logposts)
  postbank[,,iter] <- logposts
  for(j in 1:Q){
    progeny[,,j,iter] <- matlist[[j]]
  }
  posts <- pre.curr.lpsums <- colSums(logposts)
  bestN <- order(posts, decreasing=TRUE)[1:N]
  curr.lpsums <- pre.curr.lpsums[bestN]
  lpsums[iter,] <- curr.lpsums
  iter <- iter + 1
  if(iter<=iters){
    for(j in 1:length(bestN)){
      cnarray[,,j,iter] <- matlist[[bestN[j]]]
    }
  }
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
#accs[iters,]
#lpsums[iters,]

lpmaxes <- sapply(1:iters, function(j){max(lpsums[j,])})
lpmeds <- sapply(1:iters, function(j){median(lpsums[j,])})
map.indices <- sapply(1:iters, function(j){which.max(lpsums[j,])})
resids <- sapply(1:iters, function(j){sum((t(exparray[,,map.indices[j],j]) - data[,1:2])^2)})
accs.best <- sapply(1:iters, function(j){accs[j,which.max(lpsums[j,])]})
plot(lpmaxes, main='Max Log Posteriors', xlab='Iteration', ylab='Max Log Posterior')
plot(lpmeds, main='Median Log Posteriors', xlab='Iteration', ylab='Median Log Posterior')
plot(accs.best, main='Accuracy', xlab='Iterations', ylab='Overall Accuracy')
plot(resids, main='Sum of Squared Residuals' , xlab='Iteration', ylab='Sum of R^2')
plotfun1 <- function(i){
  best.acc.row <- sapply(1:iters, function(j){seg.acc.array[i,which.max(lpsums[j,]),j]})
  plot(best.acc.row, main=paste('Row ', i, ' Accuracy', sep=''), xlab='Iteration', 
       ylab='Row Accuracy')
}
sapply(1:nrow(data), plotfun1)
