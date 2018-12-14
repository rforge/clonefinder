#KRC: source('parOpt/clonefinder/functions.R')
library(CloneFinder)
library(gtools)
library(sciClone)
library(phylobase)
library(expands)
library(irr)

format.sims <- function(sim, data, kmax){
  mutids <- unique(unlist(lapply(1:length(sim$clones), function(x){sim$clones[[x]]$seq$mut.id})))
  if(is.null(mutids)){
    mutids <- numeric(0)
  }
  markers <- data$cn.data$markers
  indices <- which(sim$clones[[1]]$cn$seg %in% data$cn.data$seg)
  clones <- lapply(1:length(sim$clones),function(x){
    mutated <- rep(0, length(mutids))
    if(length(mutids)>0){
      if(length(mutids)>0){
        for(j in 1:length(mutids)){
          if(mutids[j] %in% sim$clones[[x]]$seq$mut.id){
            index <- which(sim$clones[[x]]$seq$mut.id==mutids[j])
            mutated[j] <- sim$clones[[x]]$seq$mutated.copies[index]
          }else{
            mutated[j] <- 0
          }
        }
      } 
    }
    mutated <- data.frame('mutid'=mutids, 'mutated.copies'=mutated)
    segs <- data.frame('seg'=sim$clones[[x]]$cn$seg[indices], 'markers'=markers,
              'CN'=sim$clones[[x]]$cn$A[indices] + sim$clones[[x]]$cn$B[indices])
    list('psi'=sim$psi[x], 'cn'=segs, 'mutated'=mutated)
  })
  ord <- sort(sapply(1:length(clones),function(x){clones[[x]]$psi}), index.return=TRUE, decreasing=TRUE)$ix
  clones <- clones[ord]
  if(length(clones)<kmax){
    nullclone <- list('psi'=0, 'cn'=data.frame('seg'=sim$clones[[1]]$cn$seg[indices], 'markers'=markers, 
        'CN'=rep(0, length(markers))), 'mutated'=data.frame('mutid'=mutids, 'mutated.copies'=rep(0, length(mutids))))
    clones <- c(clones, lapply(1:(kmax-length(clones)),function(x){nullclone}))
  }
  clones
}

standardize <- function(res, alg, mutdata, cndata, kmax){
  if(alg=='exp'){
    res <- res$res
    psi <- res$finalSPs[,1]
    ancestors <- res$finalSPs[,5]
    ancestors[which(ancestors<psi)] <- NA
    NAs <- which(is.na(ancestors))
    nMuts <- res$finalSPs[NAs,5]
    maxna <- which.max(NAs)
    exclude <- c(NAs[-maxna])
    exclude <- c(exclude, which(ancestors %in% psi[NAs]))
    psi <- psi[-exclude]
    ancestors <- ancestors[-exclude]
    mutids <- as.vector(res$dm[,which(colnames(res$dm)=='mutid')])
    segs <- unlist(sapply(1:length(mutids),function(k){mutdata$seg[which(mutdata$mut.id==mutids[k])]}))
    segs.unique <- unique(segs)
    markers <- cndata$markers[segs.unique]
    if(is.null(res$dm)){
      clones <- list(list('psi'=c(1),'cn'=data.frame('seg'=numeric(0),'markers'=numeric(0),'CN'=numeric(0)),
                          'mutated'=data.frame('mutid'=numeric(0),'mutated.copies'=numeric(0))))
    }else if(ncol(res$dm)>15){
      clones <- lapply(1:(length(psi)+1),function(x){
        if(x<=length(psi)){
          cnvec <- res$dm[,12]
          cnvec[which(res$dm[,15+x]==1)] <- res$dm[which(res$dm[,15+x]==1),13]
          cnvec <- sapply(1:length(segs.unique),function(y){round(mean(cnvec[which(segs==segs.unique[y])]))})
          if(length(segs.unique)==0){
            cnvec <- numeric(0)
          }
          mutcopies <- res$dm[,15+x]
          if(length(mutids)==0){
            mutcopies <- numeric(0)
          }
          output <- list('psi'=unname(psi[x]-sum(psi[which(ancestors==psi[x])])), 'cn'=data.frame('seg'=segs.unique, 'markers'=markers,                                                                     
                                                                                                  'CN'=cnvec), 'mutated'= data.frame('mutid'=mutids,'mutated.copies'=mutcopies))
        }else{
          cnvec <- sapply(1:length(segs.unique), function(y){round(mean(res$dm[which(segs==segs.unique[y]),12]))})
          if(length(segs.unique)==0){
            cnvec <- numeric(0)
          }
          output <- list('psi'=unname(1-sum(psi[which(is.na(ancestors))])), 'cn'=data.frame('seg'=segs.unique,'markers'=markers,
                                                                                            'CN'=cnvec), 
                         'mutated'=data.frame('mutid'=mutids, 'mutated.copies'=rep(0, length(mutids))))
        }
        output
      })
    }else{
      psi <- psi[which.max(sapply(1:length(psi),function(x){length(which(res$dm[,10]==psi[x]))}))]
      clones <- list(list('psi'=psi,'cn'=data.frame('seg'=segs.unique,'markers'=markers,
                                                    'CN'=sapply(segs.unique, function(x){round(mean(res$dm[which(segs==x),13]))})),
                          'mutated'=data.frame('mutid'=mutids, 'mutated.copies'=rep(1, length(mutids)))), list('psi'=unname(1-psi), 
                                                                                                               'cn'=data.frame('seg'=segs.unique,'markers'=markers,
                                                                                                                               'CN'=sapply(segs.unique, function(x){round(mean(res$dm[which(segs==x),12]))})), 
                                                                                                               'mutated'=data.frame('mutid'=mutids, 'mutated.copies'=rep(0, length(mutids)))))
    }
  }else if(alg=='sc'){
    sc <- res$res
    psi <- rep(0, kmax)
    psi.sc <- sort(sc@clust$pi, decreasing=TRUE)
    psi[1:length(psi.sc)] <- psi.sc
    merged <- sc@vafs.merged
    mutids <- merged[,6]
    segs <- unlist(sapply(1:length(mutids),function(k){mutdata$seg[which(mutdata$mut.id==mutids[k])]}))
    segs.unique <- unique(segs)
    markers <- cndata$markers[segs.unique]
    clones <- lapply(1:length(psi.sc),function(x){
      muts <- rep(0, nrow(merged))
      muts[which(merged$cluster==x)] <- 1
      mutated <- data.frame('mutid'=mutids, 'mutated.copies'=muts)
      CN <- sapply(segs.unique, function(y){
        round(mean(merged[which(segs==y),8]))
      })
      cn <- data.frame('seg'=segs.unique, 'markers'=markers, 'CN'=CN)
      list('psi'=psi[x], 'cn'=cn, 'mutated'=mutated)
    })
  }else if(alg=='cf'){
    if(!is.null(res$res)){
      res <- res$res
    }
    clones <- lapply(1:length(res$psi), function(x){
      if(res$psi[x]==0){
        cnvec <- rep(0, nrow(cndata))
      }else{
        cnvec <- rep(2, nrow(cndata))
      }
      if(class(res)=='matrix'){
        cnvec[as.numeric(rownames(res$A))] <- res$A[,x] + res$B[,x]
      }
      if(class(res$mutated)=='matrix'){
        out <- list('psi'=res$psi[x], 'cn'=data.frame('seg'=cndata$seg,'markers'=cndata$markers,'CN'=cnvec), 
                    'mutated'=data.frame('mutid'=res$filtered.data$mutdata.filt$mut.id, 'mutated.copies'=res$mutated[,x]))
      }else{
        out <- list('psi'=res$psi[x], 'cn'=data.frame('seg'=cndata$seg,'markers'=cndata$markers,'CN'=cnvec), 
                    'mutated'=data.frame('mutid'=numeric(0), 'mutated.copies'=numeric(0)))
      }
      out
    })
  }
  ord <- sort(sapply(1:length(clones),function(x){clones[[x]]$psi}), index.return=TRUE, decreasing=TRUE)$ix
  clones <- clones[ord]
  if(length(clones)<kmax){
    nullclone <- list('psi'=0, 'cn'=data.frame('seg'=segs.unique, 'markers'=markers, 'CN'=rep(0,length(markers))),
                      'mutated'=data.frame('mutid'=mutids, 'mutated.copies'=rep(0, length(mutids))))
    clones <- c(clones, lapply(1:(kmax-length(clones)),function(x){nullclone}))
  }
  psis <- sapply(1:length(clones),function(j){clones[[j]]$psi})
  for(j in which(psis==0)){
    clones[[j]]$cn$CN <- rep(0,length(clones[[j]]$cn$CN))
    clones[[j]]$mutated$mutated.copies <- rep(0,length(clones[[j]]$mutated$mutated.copies))
  }
  clones
}

concordance <- function(x, y){
  if(class(x)=='matrix'){
    x <- as.vector(x)
  }
  if(class(y)=='matrix'){
    y <- as.vector(y)
  }
  if(sd(x)==0 | sd(y)==0){
    concord <- length(which(x==y))/length(x)
  }else{
    rho <- cor(x,y)
    mu_x <- mean(x)
    mu_y <- mean(y)
    sigma_x <- sd(x)
    sigma_y <- sd(y)
    concord <- 2*rho*sigma_x*sigma_y/(sigma_x^2 + sigma_y^2 + (mu_x - mu_y)^2)
  }
  concord
}

assessAlgs <- function(resultpaths, algs, simpath, datapath){
  sim <- get(load(simpath))
  dat <- get(load(datapath))
  if(!is.null(sim$tumor)){
    sim <- sim$tumor
  }
  if(!is.null(dat$dat)){
    dat <-dat$dat
  }
  cndata <- dat$cn.data
  tr <- format.sims(sim, dat, kmax=5)
  temp <- get(load(resultpaths[[which(unlist(algs)=='cf')]]))$res
  mutdata.filt <- temp$filtered.data$mutdata.filt
  standard <- lapply(1:length(resultpaths),function(j){
    standardize(get(load(resultpaths[[j]])), algs[[j]], mutdata.filt, cndata, kmax=5)
  })
  common.seg <- standard[[1]][[1]]$cn$seg
  common.mut <- standard[[1]][[1]]$mutated$mutid
  if(length(standard)>1){
    for(j in 2:length(standard)){
      common.seg <- intersect(common.seg, standard[[j]][[1]]$cn$seg)
      common.mut <- intersect(common.mut, standard[[j]][[1]]$mutated$mutid)
    } 
  }
  absent.in.true <- common.mut[!common.mut %in% tr[[1]]$mutated$mutid]
  absent.in.res <- tr[[1]]$mutated$mutid[!tr[[1]]$mutated$mutid %in% common.mut]
  false.pos.muts <- data.frame('mutid'=absent.in.true,'mutated.copies'=rep(0,length(absent.in.true)))
  false.neg.muts <- data.frame('mutid'=absent.in.res,'mutated.copies'=rep(0,length(absent.in.res)))
  lapply(1:length(standard),function(j){
    true <- tr
    res <- standard[[j]]
    for(k in 1:5){
      res[[k]]$cn <- res[[k]]$cn[which(res[[k]]$cn$seg %in% common.seg),]
      res[[k]]$cn$CN[is.na(res[[k]]$cn$CN)] <- 2
      res[[k]]$mutated <- res[[k]]$mutated[which(res[[k]]$mutated$mutid %in% common.mut),]
      res[[k]]$mutated <- rbind(res[[k]]$mutated,false.neg.muts)
      true[[k]]$cn <- true[[k]]$cn[which(true[[k]]$cn$seg %in% common.seg),]
      true[[k]]$mutated <- rbind(true[[k]]$mutated,false.pos.muts)
      res[[k]]$cn <- res[[k]]$cn[with(res[[k]]$cn,order(res[[k]]$cn$seg)),]
      res[[k]]$mutated <- res[[k]]$mutated[with(res[[k]]$mutated,order(res[[k]]$mutated$mutid)),]
      true[[k]]$cn <- true[[k]]$cn[with(true[[k]]$cn,order(true[[k]]$cn$seg)),]
      true[[k]]$mutated <- true[[k]]$mutated[with(true[[k]]$mutated,order(true[[k]]$mutated$mutid)),]
    }
    tr.psis <- sapply(1:5,function(k){tr[[k]]$psi})
    psi.diffs <- sapply(1:5,function(k){res[[k]]$psi - true[[k]]$psi})
    cndiffs <- t(sapply(1:5,function(k){res[[k]]$cn$CN - true[[k]]$cn$CN}))
    mutdiffs <- t(sapply(1:5,function(k){res[[k]]$mutated$mutated.copies - true[[k]]$mutated$mutated.copies}))
    if(nrow(cndiffs)!=5){
      cndiffs <- t(cndiffs)
    }
    if(length(common.seg)==0){
      weighted.cndiffs <- NA
      concord.cn <- NA
    }else{
      weighted.cndiffs <- t(tr.psis*cndiffs)
      concord.cn <- concordance(sapply(1:5,function(j){true[[j]]$cn$CN}),sapply(1:5,function(j){res[[j]]$cn$CN}))
    }
    if(nrow(mutdiffs)!=5){
      mutdiffs <- t(mutdiffs)
    }
    if(length(common.mut)==0){
      weighted.mutdiffs <- NA
      concord.mut <- NA
    }else{
      weighted.mutdiffs <- t(tr.psis*mutdiffs)
      concord.mut <- concordance(sapply(1:5,function(j){true[[j]]$mutated$mutated.copies}),sapply(1:5,function(j){res[[j]]$mutated$mutated.copies}))
    }
    kdiff <- length(which(sapply(1:5,function(k){standard[[j]][[k]]$psi})>0)) - length(which(sapply(1:5,function(k){true[[k]]$psi})>0))
    concord.psi <- concordance(sapply(1:5,function(k){true[[k]]$psi}),sapply(1:5,function(k){standard[[j]][[k]]$psi}))
    list('k.dist'=(kdiff^2)^.5,'psi.dist'=sum(psi.diffs^2)^.5,'cn.dist'=sum(weighted.cndiffs^2)^.5,'mut.dist'=sum(weighted.mutdiffs^2)^.5,
         'kdiff'=kdiff,'avg.psidiff'=mean(psi.diffs),'avg.cndiff'=mean(weighted.cndiffs),'avg.mutdiff'=mean(weighted.mutdiffs),
         'concord.psi'=concord.psi,'concord.cn'=concord.cn,'concord.mut'=concord.mut,'true'=true,'res'=res)
  })
}

#####################################################################
algs <- list('cf')
simpaths <- sapply(1:300,function(j){paste('parOpt/sim/sims-mix/mixsim-',j,'.rda',sep='')})
datapaths <- sapply(1:300,function(j){paste('parOpt/dat/data-mix/mixdat-', j, '.rda', sep='')})
resultpaths <- lapply(1:300,function(j){paste('parOpt/res/res-mix/res-cf-mix/res-cf-mix-', j, '.rda', sep='')})
assessments.mix <- lapply(1:300, function(i){assessAlgs(resultpaths[[i]],algs,simpaths[i],datapaths[i])})
save(assessments.mix, file='parOpt/assess/assess-mix/assessments-mix-26.rda')

algs <- list('cf','sc','exp')
simpaths <- sapply(1:300,function(j){paste('parOpt/sim/sims-1/sim-1-',j,'.rda',sep='')})
datapaths <- sapply(1:300,function(j){paste('parOpt/dat/data-1/dat-1-', j, '.rda', sep='')})
resultpaths <- lapply(1:300,function(j){list(paste('parOpt/res/res-1/res-cf-1/res-cf-1-', j, '.rda', sep=''),
                                             paste('parOpt/res/res-1/res-sc-1/res-sc-1-', j, '.rda', sep=''),
                                             paste('parOpt/res/res-1/res-exp-1/res-exp-1-', j, '.rda', sep=''))})
assessments.1 <- lapply(1:300, function(i){assessAlgs(resultpaths[[i]],algs,simpaths[i],datapaths[i])})
save(assessments.1, file='parOpt/assess/assess-1/assessments-1-4.rda')

algs <- list('cf','sc','exp')
simpaths <- sapply(1:300,function(j){paste('parOpt/sim/sims-2/sim-2-',j,'.rda',sep='')})
datapaths <- sapply(1:300,function(j){paste('parOpt/dat/data-2/dat-2-', j, '.rda', sep='')})
resultpaths <- lapply(1:300,function(j){list(paste('parOpt/res/res-2/res-cf-2/res-cf-2-', j, '.rda', sep=''),
                                             paste('parOpt/res/res-2/res-sc-2/res-sc-2-', j, '.rda', sep=''),
                                             paste('parOpt/res/res-2/res-exp-2/res-exp-2-', j, '.rda', sep=''))})
assessments.2 <- lapply(1:300, function(j){assessAlgs(resultpaths[[j]],algs,simpaths[j],datapaths[j])})
save(assessments.2, file='parOpt/assess/assess-2/assessments-2-5.rda')

algs <- list('cf')
simpaths <- sapply(1:300,function(j){paste('parOpt/sim/sims-3/sim-3-',j,'.rda',sep='')})
datapaths <- sapply(1:300,function(j){paste('parOpt/dat/data-3/dat-3-', j, '.rda', sep='')})
resultpaths <- lapply(1:300,function(j){paste('parOpt/res/res-3/res-cf-3/res-cf-3-', j, '.rda', sep='')})
assessments.3 <- lapply(1:300, function(i){assessAlgs(resultpaths[[i]],algs,simpaths[i],datapaths[i])})
save(assessments.3, file='parOpt/assess/assess-3/assessments-3-3.rda')

