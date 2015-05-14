library(gtools)
library(combinat)
names <- as.vector(read.table("D:\\Mark\\Simulations2\\namelist.txt")[,1])
callfun <- function(i){
  string <- names[i]
  #dat <- get(load(paste("F:\\Simulations\\", string, ".rda", sep="")))
  dat <- get(load(paste("D:\\Mark\\Simulations2\\", string, ".rda", sep="")))
  list(dat[[1]], dat[[3]], dat[[4]])
}
cset <- t(combn(1:10, 3))
start <- t(sapply(1:nrow(cset), function(i){sort(cset[i,]/sum(cset[i,]))}))
x <- 1:1000

f <- function(sim){
  set <- callfun(sim)
  truepsis <- set[[1]]
  phiset <- set[[2]]
  zeds <- set[[3]]  
  lengths <- sapply(1:1000, function(i){length(which(phiset[i,]>=.05))})
  subsets <-list(which(lengths==1),which(lengths==2),which(lengths==3),which(lengths==4),which(lengths==5))
  nclones <- length(truepsis)
  nonones <- phiset[sort(unique(c(subsets[[2]], subsets[[3]], subsets[[4]], subsets[[5]]))),]
  i.non <- sort(c(subsets[[2]], subsets[[3]], subsets[[4]], subsets[[5]]))
  i.one <- subsets[[1]]
foo <- function(v){
  combos <- sort(unique(unlist(sapply(1:length(v), function(i){combn(v, i, sum)}))))
  m <- t(sapply(1:length(v), function(i){sort(c(v[i], 1-v[i]), decreasing=TRUE)}))
  rowmatch <- function(row){
    l <- which(sapply(1:length(subsets), function(i){(row %in% subsets[[i]])}))
    resfun <- function(j){
      if(l>=3){
        vec <- sort(v, decreasing=TRUE)
      }else{
        vec <- m[j,]
      }
      
      res <- sum((sort(phiset[row,], decreasing=TRUE)[1:length(vec)] - vec)^2)
      if(l>length(vec)){
        res <- res + sum((sort(phiset[row,], decreasing=TRUE)[(length(vec)+1):l])^2)
      }
      res
    }
    mat <- NULL
    if(l>=3){
      pick <- sort(v, decreasing=TRUE)
      mat <- match(sort(phiset[row,], decreasing=TRUE)[1:length(v)], phiset[row,])
    }else{
      pick <- m[which.min(sapply(1:nrow(m), resfun)),]
      mat <- match(sort(phiset[row,], decreasing=TRUE)[1:l], phiset[row,])
    }
    blank <- rep(0, 5)
    blank[mat] <- pick
    blank
  }
  trial <- t(sapply(i.non, rowmatch))
}
q <- which.min(sapply(1:nrow(start), function(k){prod((5*rowMeans((foo(start[k,]) - nonones)^2))^.01)}))
phiest <- foo(start[q,])
guess <- sort(start[q,])
combos <- unique(unlist(sapply(1:length(guess), function(i){combn(guess, i, sum)})))
sums <- unique(unlist(sapply(2:length(guess), function(i){combn(guess, i, sum)})))
onefun <- function(i){
  zeroes <- rep(0, 5)
  zeroes[which.max(phiset[i,])] <- 1
  zeroes
}
ones <- t(sapply(i.one, onefun))
clonegen <- function(v){
  psi <- v
  innerfun <- function(i){
   combo <- which(!phiest[i,] %in% guess & phiest[i,]!=0)
   i1 <- which(phiest[i,]==psi)
   if (length(i1)==0){
    i1 <- combo
   }
   zeroes <- rep(0, 5)
   zeroes[i1] <- 1
   zeroes
  }
  imat <- t(sapply(1:nrow(phiest), innerfun))
  blank <- matrix(rep(0, 5000), nrow=1000, ncol=5)
  blank[i.non,] <- imat
  blank[i.one,] <- ones
  blank
}
Zs <- lapply(guess, clonegen)
list(Zs, "accuracy"=sapply(1:length(guess), function(j){1-sum(zeds[[j]]!=Zs[[j]])/(2*1000)}), 
     "true psis"=truepsis, "estimates"=guess)
}

obj <- f(115)
#Works well for all but the 150s, for which it is only about 90% accurate; one psi = .5, wich causes problems.

resids <- sapply(1:nrow(start), function(i){sum((start[i,]-truepsis)^2)})
i.best <- which.min(resids)
truebest <- start[i.best,]
guess <- obj[[4]]

real <- foo(truebest)
picked <- foo(guess)
