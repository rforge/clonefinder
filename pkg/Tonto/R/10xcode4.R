#discrete scaled beta distribution
rsb <- function(n, alpha, beta, scalar, minim){
  betas <- rbeta(n, alpha, beta)
  scaled <- round(betas*scalar)
  scaled + minim
}

dsb <- function(x, alpha, beta, scalar, minim, log=FALSE){
  betas <- (x - minim)/scalar
  dbeta(betas, alpha, beta, log=log)
}

#We'll assume a uniform distribution for now:
LRs <- seq(from=minim, to=maxim, by=2)
probs <- dunif(LRs, minim, maxim)
probs <- probs/sum(probs)
psiset <- seq(from=0, to=1, length=101)

#N is the number of lengths to sample from the presumed molecule length distribution.
logLike <- function(data, psi, dist='unif', pars, cn, seglen, N){
  LA=cn*seglen
  if(dist=='unif'){
    minim <- pars$minim
    maxim <- pars$maxim
    midpoint <- minim + .5*(maxim-minim)
    L <- midpoint
    if(L >= LA){
      kset <- c(1:(L - LA), rep((L - LA), (L - LA - 1)), (L - LA):1)
    }else{
      kset <- c(1:L, rep(L, LA - L - 1), L:1)
    }
    pset <- sapply(1:nrow(data), function(i){
      n <- data[i,1]
      p <- (kset/L)^n
      (mean(p)*(1-psi) + (1*psi))
    })
  }else if(dist=='sn'){
    alpha <- pars$alpha
    omega <- pars$omega
    tau <- pars$tau
    x <- pars$x
    Ls <- round(rsn(20, xi=x, omega=omega, alpha=alpha, tau=tau))
    probs <- dsn(Ls, x, omega, alpha, tau)
    probs <- probs/sum(probs)
    pset <- sapply(1:nrow(data), function(i){
      sum(sapply(1:length(Ls), function(j){
        L <- Ls[j]
        if(L >= LA){
          kset <- c(1:(L - LA), rep((L - LA), (L - LA)), (L - LA):1)
        }else{
          kset <- c(1:L, rep(L, LA - L), L:1)
        }
        n <- data[i,1]
        p <- (kset/L)^n
        #(dbinom(n, n, prob=p)*(1 - psi) + psi)*probs[j]
        (mean(p)*(1-psi) + (1*psi))*probs[j]
      }))
    })
  }else if(dist=='beta'){
    alpha <- pars$alpha
    beta <- pars$beta
    scalar <- pars$scalar
    minim <- pars$minim
    Ls <- rsb(N, alpha, beta, scalar, minim)
    probs <- dsb(Ls, alpha, beta, scalar, minim)
    probs <- probs/sum(probs)
    pset <- sapply(1:nrow(data), function(i){
      sum(sapply(1:length(Ls), function(j){
        L <- Ls[j]
        if(L >= LA){
          kset <- c(1:(L - LA), rep((L - LA), (L - LA)), (L - LA):1)
        }else{
          kset <- c(1:L, rep(L, LA - L), L:1)
        }
        n <- data[i,1]
        p <- (kset/L)^n
        #(dbinom(n, n, prob=p)*(1 - psi) + psi)*probs[j]
        (mean(p)*(1-psi) + (1*psi))*probs[j]
      }))
    })
  }
  mu <- sum(pset)
  obs <- length(which(data[,1]==data[,2]))
  c(dpois(obs, lambda=mu, log=TRUE), psi, mu, obs, psi*nrow(data))
}

###Actual data:
GS1015 <- read.csv('C:\\Users\\Mark\\OneDrive - The Ohio State University\\Data_10x\\GS1015.csv')
GS1174 <- read.csv('C:\\Users\\Mark\\OneDrive - The Ohio State University\\Data_10x\\GS1174.csv')
GS1194 <- read.csv('C:\\Users\\Mark\\OneDrive - The Ohio State University\\Data_10x\\GS1194.csv')

#Estimating length distribution:
totals <- c(GS1015$total, GS1174$total, GS1194$total)
est <- totals/.7
#plot(density(est), xlim=c(9, 200))
sd.guess <- sd(est)

#Not sure; instead, let's just set mu = 65 and assume the distribution skews to the right.
dat1015 <- data.frame('size'=GS1015$total, K=GS1015$total - GS1015$A - GS1015$C, 'V'=GS1015$V)
dat1015 <- dat1015[which(dat1015$V>0),1:2]
dat1174 <- data.frame('size'=GS1174$total, K=GS1174$total - GS1174$A - GS1174$C, 'V'=GS1174$V)
dat1174 <- dat1174[which(dat1174$V>0),1:2]
dat1194 <- data.frame('size'=GS1194$total, K=GS1194$total - GS1194$A - GS1194$C, 'V'=GS1194$V)
dat1194 <- dat1194[which(dat1194$V>0),1:2]

datasets.real <- list(dat1015, dat1174, dat1194)
nameset <- c('GS1015', 'GS1174', 'GS1194')
####End of stuff to load in...


###Probability distribution for long molecule length, assuming a scaled beta distribution:
alpha.range <- seq(from=2, to=20, length=6)
beta.range <- seq(from=1, to=41, length=6)
scalar.range <- seq(from=80, to=150, length=6)
minim <- 20
partab <- expand.grid(alpha.range, beta.range, scalar.range)
partab <- cbind(partab, rep(minim, nrow(partab)))
parnames <- c('alpha', 'beta', 'scalar', 'minim')
partab <- as.data.frame(partab)
colnames(partab) <- parnames
mu.m <- 40
mu.s <- 5
sigma.m <- 1
sigma.s <- .5
parprior <- sapply(1:nrow(partab), function(i){
  mu <- partab$minim[i] + partab$scalar[i]*partab$alpha[i]/(partab$alpha[i] + partab$beta[i])
  sigma <- partab$scalar[i]*((partab$alpha[i]*partab$beta[i])/(((partab$alpha[i] + 
                                                                   partab$beta[i] + 1))*(partab$alpha[i] + partab$beta[i])^2))^.5
  dnorm(mu, mu.m, sigma.m)*dnorm(sigma, mu.s, sigma.s)
})
parprior <- parprior/sum(parprior)
#plot(parprior, type='l')
partab$prior <- parprior

lrange <- 20:120
pre <- 1:(lrange[1]-1)
len.probs <- sapply(lrange, function(i){
  sum(sapply(1:nrow(partab),function(j){
    dsb(i, partab$alpha[j], partab$beta[j], partab$scalar, partab$minim, log=FALSE)*partab$prior[j]
  }))
})
plot(c(pre, lrange), c(rep(0, length(pre)),len.probs/sum(len.probs)), xlab='Length', ylab='Probability', 
     main='Length Distribution', type='l')

#Compute posterior probability of psi (episomal fraction)
compPost <- function(k, q, dist='beta'){
  psi <- psiset[k]
  t1 <- Sys.time()
  res <- lapply(1:nrow(partab), function(j){
    pars <- as.list(partab[j,])
    names(pars) <- parnames
    llres <- logLike(data=datasets.real[[q]], psi=psi, dist=dist, pars=pars, cn=5, seglen=50, N=20)
    list('psi.prob'=llres[1]*partab$prior[j])
  })
  t2 <- Sys.time()
  tdiff <- t2 - t1
  lpp <- unlist(lapply(1:length(res), function(j){res[[j]][[1]]}))
  pp <- exp(lpp)
  ppsum <- sum(pp)
  list('ppsum'=ppsum)
}

#q is the sample index; partab is the table of length distribution parmaeters; parnames is the column names;
#int is the interval for the CI (default 95% CI); path is the location to save the image.
compCI <- function(q, psiset, partab, parnames, int=.95, path){
  res <- lapply(1:length(psiset), compPost, q=q, dist='beta')
  posts <- unlist(sapply(1:length(res), function(i){res[[i]][[1]]}))
  posts <- posts/sum(posts)
  maxindex <- which.max(posts)
  maxpost <- max(posts)
  remainder <- (int - maxpost)/2
  right.sum <- sum(posts[(maxindex+1):length(posts)])
  left.sum <- sum(posts[(maxindex-1):1])
  if(left.sum < remainder & right.sum < remainder){
    cred.int <- maxindex
  }else{
    if(left.sum < remainder){
      left <- 1
      remainder <- int - maxpost - left.sum
      right <- maxindex + min(which(sapply((maxindex+1):length(posts), function(i){sum(posts[(maxindex+1):i])})>remainder))
    }else if(right.sum < remainder){
      right <- 101
      remainder <- int - maxpost - right.sum
      left <- maxindex - min(which(sapply((maxindex-1):1, function(i){sum(posts[(maxindex-1):i])})>remainder))
    }else{
      right <- maxindex + min(which(sapply((maxindex+1):length(posts), function(i){sum(posts[(maxindex+1):i])})>remainder))
      left <- maxindex - min(which(sapply((maxindex-1):1, function(i){sum(posts[(maxindex-1):i])})>remainder))
    }
    cred.int <- c(left, right)
  }
  png(paste(path, '\\confint-',nameset[q], '.png', sep=''))
  plot(psiset, posts, xlab='Episomal Fraction', ylab='Posterior', main='Post. Distribution with 95% CI', type='l')
  abline(v=psiset[maxindex], col='red')
  abline(v=psiset[c(cred.int)], col='darkorange')
  #legend()
  dev.off()
  list('psis'=psiset, 'posts'=posts, 'CI'=psiset[cred.int], 'interval'=int)
}

#Compute confindence intervals
mypath <- 'C:\\Users\\Mark\\OneDrive - The Ohio State University\\Images\\10x\\confints'
confint1 <- compCI(1, psiset, partab, parnames, path=mypath)
confinft2 <- compCI(2, psiset, partab, parnames, path=mypath)
confint3 <- compCI(3, psiset, partab, parnames, path=mypath)
