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
  # load the simulated data set
  set <- callfun(sim)
  truepsis <- set[[1]] # omniscient
  phiset <- set[[2]]   # previously estimated phi from the observed data
  zeds <- set[[3]]     # omniscient
  lengths <- sapply(1:1000, function(i){length(which(phiset[i,]>=.05))}) # putative number of clones per segment
  subsets <-list(which(lengths==1),
                 which(lengths==2),
                 which(lengths==3),
                 which(lengths==4),
                 which(lengths==5)) # why this data structure/algorithm?
  nclones <- length(truepsis) # omniscient
  nonones <- phiset[sort(unique(c(subsets[[2]], subsets[[3]], subsets[[4]], subsets[[5]]))),] # informative
  i.non <- sort(c(subsets[[2]], subsets[[3]], subsets[[4]], subsets[[5]]))                # non-informative
  i.one <- subsets[[1]] # index of non-informative segments
  foo <- function(v){
    combos <- sort(unique(unlist(sapply(1:length(v), function(i){combn(v, i, sum)})))) # expand psis to sums
    m <- t(sapply(1:length(v), function(i){sort(c(v[i], 1-v[i]), decreasing=TRUE)}))   # what is this?
    rowmatch <- function(row){
      # what is thsi fiunction supposed to do? Inputs? outputs? goal?

      # next line looks like a complicated way to compute lengths[row]. am I missing something?
      l <- which(sapply(1:length(subsets), function(i){(row %in% subsets[[i]])}))
      resfun <- function(j){
        # note: this function only gets called when l < 3. Why another if-else?
        if(l >= 3){ # it's probably bad practice to keep relying on external vraiables instead of passing them in as arguments
          vec <- sort(v, decreasing=TRUE)
        }else{ # what is vec for? and whow does m come into it?
          vec <- m[j,]
        }
        
        res <- sum((sort(phiset[row,], decreasing=TRUE)[1:length(vec)] - vec)^2) # sum of squared errors?
        if(l > length(vec)){
          res <- res + sum((sort(phiset[row,], decreasing=TRUE)[(length(vec)+1):l])^2) # ???
        }
        res
      }
      mat <- NULL
      if(l >= 3){
        pick <- sort(v, decreasing=TRUE)
        mat <- match(sort(phiset[row,], decreasing=TRUE)[1:length(v)], phiset[row,]) # what are we matchinjg?
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
  q <- which.min(sapply(1:nrow(start), function(k){
                          prod((5*rowMeans((foo(start[k,]) - nonones)^2))^.01)
                        }))
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
