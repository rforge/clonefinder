dbeta2 <- function(x, mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  dbeta(x, abs(alpha), abs(beta))
}

rbeta2 <- function(n, mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  rbeta(n, abs(alpha), abs(beta))
}

dgamma2 <- function(x, mode, sigma, log=FALSE){
  mu <- .5*(mode + (mode^2 + 4*sigma^2)^.5)
  scale <- (sigma^2)/mu
  shape <- mu/scale
  dgamma(x, shape=shape, scale=scale, log=log)
}

rgamma2 <- function(n, mode, sigma, log=FALSE){
  mu <- .5*(mode + (mode^2 + 4*sigma^2)^.5)
  scale <- (sigma^2)/mu
  shape <- mu/scale
  log(rgamma(n, shape=shape, scale=scale))
}

logprior.j <- function(psi){
  log(ddirichlet(psi, rep(1/2, length(psi))))
}

logprior.cn1 <- function(A, B, prob){
  sum(dgeom(A, prob, log=TRUE)) + sum(dgeom(B, prob, log=TRUE))
}

logprior.cn2 <- function(A, B, theta){
  sum(dexp(A, theta, log=TRUE)) + sum(dexp(B, theta, log=TRUE))
}

logprior.k <- function(psi, theta2){
  dexp(length(which(psi!=0)), theta2, log=TRUE)
}

ll.gamma <- function(z, dat, psivec, marks, sigma0, nclone){
  a <- z[1:nclone]
  b <- z[(nclone+1):(2*nclone)]
  sigma <- sigma0/(marks^.5)
  mode1 <- sum(a*psivec)
  mode2 <- sum(b*psivec)
  dgamma2(unlist(dat[1]), mode1, sigma, log=TRUE) + dgamma2(unlist(dat[2]), mode2, sigma, log=TRUE)
}

lp <- function(dat, marks, z, psivec, theta, sigma0, nclone){
  term1 <- logprior.j(psivec)
  A <- z[1:nclone]
  B <- z[(nclone+1):(2*nclone)]
  term2 <- logprior.cn2(A, B, theta) + ll.gamma(z, dat, psivec, marks, sigma0, nclone=nclone)
  term1 + term2
}

likely.seg <- function(datavec, marks, ABvec, psi, sigma0){
  datavec <- unlist(datavec)
  term1 <- dnorm(sum(ABvec[1:(length(ABvec)/2)]*psi) - datavec[1], 0, sigma0/(marks^.5), log=TRUE)
  term2 <- dnorm(sum(ABvec[(length(ABvec)/2+1):(length(ABvec))]*psi) - datavec[2], 0, 
                 sigma0/(marks^.5), log=TRUE)
  term1 + term2
}

prior.seg <- function(vec, psi, theta1, theta2){
  term1 <- log(ddirichlet(psi, rep(1/2, length(psi))))
  term2 <- sum(dgeom(abs(vec - 1), theta1, log=TRUE))
  term3 <- dgeom(length(psi), theta2, log=TRUE)
  term1 + term2 + term3
}

findAB <- function(data, markers, psi, N){
  mat <- matrix(sample(0:5, 2*length(psi)*N, replace=TRUE), nrow=N, ncol=2*length(psi))
  map <- matrix(rep(NA, nrow(data)*ncol(mat)), nrow=nrow(data), ncol=ncol(mat))
  for(i in 1:nrow(data)){
    lls <- unname(sapply(1:nrow(mat), function(j){likely.seg(data[i,], markers[i], mat[j,], 
                                                             psi, sigma0=.1)}))
    imax <- which.max(lls)
    cns <- sapply(1:ncol(mat), function(k){
      mns <- sapply(0:5, function(l){mean(lls[which(mat[,k]==l)])})
      which.max(mns) - 1
    })
    map[i,] <- cns
  }
}

findAB2 <- function(data, markers, psi){
  map <- matrix(rep(NA, nrow(data)*ncol(mat)), nrow=nrow(data), ncol=ncol(mat))
  maxposts <- rep(NA, nrow(data))
  postsets <- list(NULL)
  for(i in 1:nrow(data)){
    lls <- unname(sapply(1:nrow(mat), function(j){likely.seg(data[i,], markers[i], mat[j,], 
                                                             psi, sigma0=1)}))
    lps <- sapply(1:nrow(mat), function(j){prior.seg(unlist(mat[j,]), psi, theta1=.5, 
                                                     theta2=.5)})
    lpost <- lls + lps
    #plot(lpost)
    #abline(v=truebest, col='red')
    imax <- which.max(lpost)
    #cns <- sapply(1:ncol(mat), function(k){
    #  mns <- sapply(0:5, function(l){mean(lls[which(mat[,k]==l)])})
    #  which.max(mns) - 1
    #})
    cns <- unlist(mat[imax,])
    map[i,] <- cns
    #colnames(map) <- colnames(truemat)
    postsets <- c(postsets, list(lpost))
    maxposts[i] <- imax
  }
  list(map, lpost, maxposts)
}

likely.mat <- function(data, psi, mat, sigma0){
  markers <- data$markers
  expected <- data.frame('x'=sapply(1:nrow(data), function(i){sum(unlist(mat[i,1:(ncol(mat)/2)])*psi
    )}),'y'=sapply(1:nrow(data), function(i){sum(unlist(mat[i,((ncol(mat)/2)+1):ncol(mat)])*psi)}))
  dif <- data[,1:2] - expected
  sum(sapply(1:nrow(dif), function(i){dnorm(unlist(dif[i,]), 0, sigma0/markers, log=TRUE)}))
}

prior.mat <- function(mat, theta){
  vec <- as.vector(mat)
  sum(dgeom(abs(vec - 1), theta, log=TRUE))
}

#######
computeSSE <- function(a, psi0, segdata) {
  err <- sum(a*psi0) - segdata
  err^2
}

posteriorAvector <- function(a, psi, sig, segdata, prior) {
  sse <- computeSSE(a, psi, segdata)
  like <- dnorm(sse, 0, sig)
  post <- like*prod(prior[1+a])
  c(a, sse=sse, like=like, post=post)
}



#######
findPsi <- function(data, AB){
  AB <- as.matrix(AB)
  M <- rbind(AB[,1:(ncol(AB)/2)], AB[,(ncol(AB)/2+1):ncol(AB)])
  dat <- as.matrix(c(data[,1], data[,2]))
  psi <- as.vector(solve(t(M) %*% M) %*% t(M) %*% dat)
  psi <- psi/sum(psi)
  psi
}
#(M(T)M)^(-1)*M(T)*W; W = c(X, Y) 

bestpick <- function(z, psi, markers, theta, sigma0, N=1000){
  prior <- dgeom(abs(1-0:7), prob=theta)
  postmat <- sapply(1:N, function(q){
    vec <- sample(0:7, replace=TRUE, prob=prior)
    sigma <- sigma0/markers
    post <- dnorm(vec*psi - z, 0,  sigma, log=TRUE) + sum(log(prior[c(vec+1)]))
    c(vec, post)
  })
  newposts <- t(sapply(1:length(psi), function(j){sapply(0:7, function(x){
    if(length(which(postmat[,j]==x))>0){
      out <- sum(postmat[,ncol(postmat)][which(postmat[,j]==x)])
    }else{
      out <- 0
    }
  })}))
  sapply(1:ncol(newposts), function(j){which.max(newposts[,j])-1})
}

findMat <- function(data, markers, psilist, sigma0=1, n=100, theta=.3){
  psi <- unlist(psilist)
  mat <- sapply(1:nrow(data), function(j){
    dat <- data[j,]
    marks <- markers[j]
    x <- dat$X
    y <- dat$Y
    best <- c(bestpick(x, psi, marks, theta=theta, sigma0), bestpick(y, marks, psi, theta=theta, sigma0))
    best
  })
  mat
}

matgen <- function(cols, rows, theta, max){
  N <- cols*rows
  integers <- 0:max
  probs <- dgeom(abs(integers - 2), prob=theta)
  vec <- sample(integers, N, replace=TRUE, prob=probs)
  matrix(vec, nrow=rows, ncol=cols)
}
