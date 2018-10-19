#load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\clonefinder\\cnmodels.rda')
#load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\clonefinder\\psis.10.rda')
#load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\clonefinder\\psis.20.rda')
#load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\clonefinder\\psis.50.rda')
#load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\clonefinder\\psis.100.rda')
#load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\clonefinder\\etamat.10.rda')
#load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\clonefinder\\etamat.20.rda')
#load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\clonefinder\\etamat.50.rda')

psiOptim <- function(cndata.filt, mutdata.filt, psis, cnmodels, pars, cnmax=5, kPriors=NULL){
  if(is.null(kPriors)){
    kPriors <- dgeom((1:5)-1,prob=pars$ktheta, log=TRUE)
  }
  kPriorset <- sapply(1:nrow(psis),function(k){kPriors[length(which(psis[k,]>pars$thresh))]})
  res <- lapply(1:nrow(psis),function(k){    nonzero <- which(psis[k,]>0)
    models <- cnmodels[,nonzero]
    if(length(nonzero)==1){
      models <- as.matrix(models)
    }
    models <- as.matrix(unique(as.data.frame(models[,nonzero])))
    cnPriors <- sapply(1:nrow(models),function(j){sum(dgeom(abs(1-models[j,]), prob=pars$theta,log=TRUE))})
    index.upper <- which(sapply(1:nrow(models),function(m){length(which(models[m,]<1))==0}))
    index.lower <- which(sapply(1:nrow(models),function(m){length(which(models[m,]>1))==0}))
    mods.upper <- as.matrix(models[index.upper,])
    mods.lower <- as.matrix(models[index.lower,])
    if(length(index.upper)==1){
      mods.upper <- t(mods.upper)
    }
    if(length(index.lower)==1){
      mods.lower <- t(mods.lower)
    }
    upper.priors <- cnPriors[index.upper]
    lower.priors <- cnPriors[index.lower]
    etas.upper <- as.vector(psis[k,nonzero] %*% t(mods.upper))
    etas.lower <- as.vector(psis[k,nonzero] %*% t(mods.lower))
    etaList <- list(etas.lower,etas.upper)
    priorList <- list(lower.priors,upper.priors)
    modList <- list(mods.lower, mods.upper)
    x.numbers <- y.numbers <- rep(1,nrow(cndata.filt))
    x.numbers[which(cndata.filt$X>1)] <- 2
    y.numbers[which(cndata.filt$Y>1)] <- 2
    xmat <- t(sapply(1:nrow(cndata.filt),function(l){
      post <- dnorm(cndata.filt$X[l]-etaList[[x.numbers[l]]],0,pars$sigma0/(cndata.filt$markers[l]^.5),log=TRUE) + 
        priorList[[x.numbers[l]]]
      avec <- c(modList[[x.numbers[l]]][which.max(post),],rep(0,cnmax-length(nonzero)))
      c(avec,etaList[[x.numbers[l]]][which.max(post)],max(post))
    }))
    ymat <- t(sapply(1:nrow(cndata.filt),function(l){
      post <- dnorm(cndata.filt$Y[l]-etaList[[y.numbers[l]]],0,pars$sigma0/(cndata.filt$markers[l]^.5),log=TRUE) + 
        priorList[[y.numbers[l]]]
      bvec <- c(modList[[y.numbers[l]]][which.max(post),],rep(0,cnmax-length(nonzero)))
      c(bvec,etaList[[y.numbers[l]]][which.max(post)],max(post))
    }))
    cnRes <- list('A'=xmat[,1:cnmax],'B'=ymat[,1:cnmax],'etaA'=xmat[,cnmax+1],'etaB'=ymat[,cnmax+1],
                  'post'=sum(xmat[,cnmax+2])+sum(ymat[,cnmax+2]))
    psiPost <- cnRes$post
    if(nrow(mutdata.filt)>0){
      mutmat <- t(sapply(1:nrow(mutdata.filt),function(l){
        cn.index <- which(cndata.filt$seg==mutdata.filt$seg[l])
        if(length(cn.index)>0){
#was:          A <- A.array[cn.index,]
#was:          B <- B.array[cn.index,]
          A <- cnRes$A[cn.index,]
          B <- cnRes$B[cn.index,]
        }else{
          A <- B <- rep(1,K)
        }
        mut.models <- as.matrix(expand.grid(lapply(1:K,function(m){0:max(c(A[m],B[m]))})))
        mll <- sapply(1:nrow(mut.models), function(m){
          if(sum((A+B)*psis[k,])>0){
            p <- sum(mut.models[m,]*psis[k,])/sum(psis[k,]*(A+B))
            if(p<.0001){
              p <- .0001
            }
          }else{
            p <- .0001
          }
#Was:          log(sum(dbinom(mutdata.filt$varCounts[l],mutdata.filt$totalCounts[l],pvec)*dbeta2(pvec,p,pars$sigmaM)))
          log(sum(dbinom(mutdata.filt$varCounts[l], mutdata.filt$totalCounts[l], p)*dbeta2(p, p, pars$sigmaM)))
        })
        mprior <- sapply(1:nrow(mut.models),function(l){dgeom(length(unique(mut.models[l,which(psis[k,]>pars$thresh)]))-1,prob=pars$mtheta, log=TRUE)})
        mposts <- mll + mprior
        maxind <- which.max(mposts)
        etaM <- sum(psis[k,]*mut.models[maxind,])/sum(psis[k,]*(A+B))
        if(is.na(etaM)){
          etaM <- 0
        }
        unname(c(mut.models[which.max(mposts),],mposts[maxind],etaM))
      }))
      mutPosts <- sum(mutmat[,ncol(mutmat)-1])
      psiPost <- psiPost + sum(mutPosts)
      etaM <- mutmat[,ncol(mutmat)]
      mutated <- mutmat[,1:(ncol(mutmat)-2)]
      rownames(mutated) <- rownames(mutdata.filt)
    }else{
      mutated <- etaM <- NA
    }
    list('A'=cnRes$A,'B'=cnRes$B,'mutated'=mutated,'psiPost'=psiPost,'etaA'=cnRes$etaA,'etaB'=cnRes$etaB,'etaM'=etaM)
  })
  psiPosts <- sapply(1:length(res),function(j){res[[j]]$psiPost})
  psiPosts <- psiPosts + kPriorset
  pick <- which.max(psiPosts)
  A <- res[[pick]]$A
  B <- res[[pick]]$B
  rownames(A) <- rownames(B) <- rownames(cndata.filt)
  list('A'=A,'B'=B,'mutated'=res[[pick]]$mutated,'etaA'=res[[pick]]$etaA,'etaB'=res[[pick]]$etaB,'etaM'=res[[pick]]$etaM,'psiPosts'=psiPosts)
}

###Other functions
psiFilter <- function(psis, psi.center, range){
  diff <- t(t(psis) - psi.center)
  bin <- abs(diff) <= range
  pass <- which(sapply(1:nrow(bin),function(k){length(which(bin[k,]))==ncol(bin)}))
  list('psis'=psis[pass,], 'indices'=pass)
}

#'snp data' means cgh data unless there is no cgh data, in which case we substitute in snp info from seq data
runAlg <- function(cndata, vardata, cnmodels, psiset, pars, imputedCN=NULL){
  if(is.null(cndata)){
    if(is.null(imputedCN)){
      seqsnps <- vardata[vardata$status=='germline',]
      snpdata <- seqSeg(seqsnps, len=100, thresh=.25)
    }else{
      snpdata <- imputedCN
    }
  }else{
    snpdata <- cndata
    #denx <- density(na.omit(snpdata$X))
    #deny <- density(na.omit(snpdata$Y))
    #xrange <- which(denx$x > .7 & denx$x < 1.3)
    #yrange <- which(deny$x > .7 & deny$x < 1.3)
    #shift.X <- 1 - denx$x[xrange][which.max(denx$y[xrange])]
    #shift.Y <- 1 - deny$x[yrange][which.max(deny$y[yrange])]
    #snpdata$X <- snpdata$X + shift.X
    #snpdata$Y <- snpdata$Y + shift.Y
    #snpdata$X[snpdata$X<0] <- 0
    #snpdata$Y[snpdata$Y<0] <- 0
  }
  if(is.null(vardata)){
    seqdata <- data.frame("chr"=numeric(),"seg"=numeric(),"mut.id"=numeric(),"refCounts"=numeric(),
                          "varCounts"=numeric(),"VAF"=numeric(),"totalCounts"=numeric(),'status'=character())
    mu <- NA
  }else{
    seqdata <- vardata
    if(nrow(vardata)>0){
      read.den <- density(vardata$refCounts)
      peak <- read.den$x[which.max(read.den$y)]
      mu <- peak
      pars$mu <- mu
      pars$sigma.counts <- sd(vardata$refCounts)
    }
  }
#was:  cn.filt <- filter(snpdata, threshold=pars$thresh, cutoff=pars$cutoff)
  cn.filt <- filter(snpdata, threshold=pars$thresh)
  cndata.filt <- cn.filt$mat
  indices.cn <- cn.filt$indices
  mutdata <- seqdata[seqdata$status=='somatic',]
  mut.filt <- filter.mut(mutdata, mu=mu, threshold=3)
  mutdata.filt <- mut.filt$mat
  mutids.filt<- mut.filt$ids
  kPriors <- dgeom((1:5)-1,prob=pars$ktheta, log=TRUE)
  psis <- psiset
  kmax <- ncol(psis)
  temp <- psiOptim(cndata.filt, mutdata.filt, psis, cnmodels, pars, cnmax=5, kPriors)
  A <- temp$A
  B <- temp$B
  mutated <- temp$mutated
  etaA <- temp$etaA
  etaB <- temp$etaB
  etaM <- temp$etaM
  posts <- temp$psiPosts
  maxPost <- max(posts)
  iters <- pars$iters
  Q <- pars$Q
  iter <- 1
  logPosts <-  rep(NA, (Q)*(iters-1)+nrow(psis))
  psibank <- matrix(NA, ncol=ncol(psis),nrow=length(logPosts))
  while(iter<=iters){
    if(iter==1){
      bank.indices <- 1:nrow(psis)
    }else{
      bank.indices <- min(which(is.na(psibank[,1]))):(min(which(is.na(psibank[,1])))+Q-1)
    }
    if(max(posts)>maxPost){
      A <- temp$A
      B <- temp$B
      mutated <- temp$mutated
      etaA <- temp$etaA
      etaB <- temp$etaB
      etaM <- temp$etaM
      maxPost <- max(posts)
    }
    logPosts[bank.indices] <- posts
    psibank[bank.indices,] <- psis
    max.index <- max(bank.indices)
    postset <- logPosts[1:max.index]
    probs <- exp(postset)
    probs[which(is.infinite(probs))] <- 1
    probs <- probs/sum(probs)
    if(length(which(is.na(probs)))==length(probs)){
      min.post <- min(postset)
      adjusted <- postset - min.post
      adjusted <- adjusted + min(adjusted)
      probs <- adjusted/sum(adjusted)
      if(length(which(is.na(probs)))==length(probs)){
        probs <- rep(1/length(probs),length(probs))
      }
    }
    if(length(which(is.infinite(postset)))==length(postset)){
      probs <- rep(1/length(probs), length(probs))
    }
    indices <- sample(1:max.index, Q, replace=TRUE, prob=probs)
    alphas <- psibank[indices,]
    psis <- t(sapply(1:nrow(alphas),function(j){
      temp <- as.vector(rdirichlet(1,alpha=alphas[j,]))
      if(length(temp)<kmax){
        temp <- c(temp, rep(0, kmax-length(temp)))
      }
      sort(temp, decreasing=TRUE)
    }))
    psis <- t(sapply(1:nrow(psis),function(x){
      foo <- psis[x,]
      foo[foo<pars$thresh] <- 0
      foo/sum(foo)
    }))
    temp <- psiOptim(cndata.filt, mutdata.filt, psis, cnmodels, pars, cnmax=5, kPriors)
    posts <- temp$psiPosts
    iter <- iter + 1
  }
  max.index <- which.max(logPosts)
  psi <- psibank[max.index,]
  psi[which(psi<pars$thresh)] <- 0
  psi <- psi/sum(psi)
  res <- list('indices'=list('indices.cn' = indices.cn, 'mutids.filt' = mutids.filt),
      'data'=list('seqdata'=seqdata, 'snpdata'=snpdata),'filtered.data'=list('mutdata.filt'=mutdata.filt,'cndata.filt'=cndata.filt),
      'psibank' = psibank,'etaA'=etaA,'etaB'=etaB,'etaM'=etaM,'psi'=psi,'A'=A,'B'=B,'mutated'=mutated,'psiPosts'=logPosts,
      'pars'=pars,'max.index'=max.index)
  res
}
