tumorGen <- function(psi, rounds, nu=100, pcnv=.5, norm.contam=FALSE, cnmax=4){
  K <- length(which(psi>0))
  total.segs <- round(runif(1, 250, 500))
  segsperchr <- as.vector(rdirichlet(1, chlens/1000000))
  segsperchr <- round(segsperchr*total.segs)
  segsperchr[segsperchr<1] <- 1
  chr <- unlist(lapply(1:24, function(i){rep(i, segsperchr[i])}))
  ends <- lapply(1:24, function(i){
    lens <- c(as.vector(rdirichlet(1, rep(1, segsperchr[i]))))
    sapply(1:length(lens),function(j){sum(lens[1:j])})
  })
  ends <- lapply(1:length(ends),function(i){round(chlens[chr[i]]*ends[[i]])})
  starts <- lapply(1:length(ends),function(i){c(1,ends[[i]][1:(length(ends[[i]])-1)]+1)})
  ends <- unlist(ends)
  starts <- unlist(starts)
  cnmat <- cbind(chr, unlist(starts), unlist(ends), rep(1, length(starts)),
                 rep(1, length(starts)), 1:length(starts),rep(NA, length(starts)))
  colnames(cnmat) <- c('chr', 'start', 'end', 'A', 'B', 'seg', 'parent.index')
  startclone <- list('cn'=cnmat, 'seq'=NA)
  id.end <- 0
  clones <- list(startclone)
  while(length(clones) <= rounds) {
    parent.index <- sample(1:length(clones),1)
    parent <- clones[[parent.index]]
    if(nu>0){
      cnv <- sample(c(TRUE, FALSE), 1, prob=c(pcnv, 1-pcnv))
    }else{
      cnv <- TRUE
    }
    newclone.cn <- parent$cn
    newclone.cn[,which(colnames(newclone.cn)=='parent.index')] <- parent.index
    newclone.seq <- parent$seq
    if(cnv){
      allele <- sample(c('A', 'B'), 1)
      pool <- which(parent$cn[,which(colnames(cnmat)==allele)]>0 & parent$cn[,which(colnames(cnmat)==allele)]<cnmax)
      tochange <- sample(pool, 1)
      other.clones <- sapply(1:length(clones),function(j){clones[[j]]$cn[tochange,which(colnames(cnmat)==allele)]})
      if(max(other.clones)>1){
        delta <- 1
      }else if(min(other.clones)<1){
        delta <- -1
      }else{
        delta <- sample(c(-1,1),1)
      }
      if(parent.index>1 & length(which(!is.na(unlist(newclone.seq))))>0){
        muts.tochange <- intersect(which(newclone.seq[,i.seg]==tochange), which(newclone.seq[,i.allele]==allele))
        if(length(muts.tochange)>0){
          newclone.seq[muts.tochange,i.cn] <- as.numeric(newclone.seq[muts.tochange,i.cn]) + delta
        }
      }
      newclone.cn[tochange,which(colnames(cnmat)==allele)] <- as.numeric(newclone.cn[tochange,which(colnames(cnmat)==allele)]) +
        delta
    }
    nmuts <- round(runif(1, nu/2, 1.5*nu))
    if(nmuts>0){
      mutsperchr <- as.vector(rmultinom(1, nmuts, chlens/1000000))
      starts <- lapply(1:length(mutsperchr),function(i){unique(round(sort(runif(mutsperchr[i], 1, chlens[i]))))})
      starts <- unlist(starts)
      mut.chrs <- unlist(lapply(1:length(mutsperchr),function(i){rep(i, mutsperchr[i])}))
      mut.segs <- unlist(sapply(1:length(starts), function(i){
        col.index <- which(colnames(cnmat)=='chr')
        tab <- cnmat[which(cnmat[,col.index]==mut.chrs[i]),]
        bin1 <- which(tab[,which(colnames(tab)=='start')]<=starts[[i]])
        index <- bin1[which(tab[bin1,which(colnames(tab)=='end')]>=starts[[i]])]
        seg <- unname(tab[index,which(colnames(tab)=='seg')])
        seg
      }))
      alleles <- sample(c('A','B'), length(starts), replace=TRUE)
      mut.ids <- (id.end+1):(id.end+length(starts))
      seqmat <- cbind(mut.chrs, starts, mut.segs, mut.ids, rep(1, length(starts)), alleles)
      id.end <- max(mut.ids)
      colnames(seqmat) <- c('chr', 'start', 'seg', 'mut.id', 'mutated.copies', 'allele')
      newclone.seq <- rbind(newclone.seq, seqmat)
      rownames(newclone.seq) <- NULL
    }else{
      newclone.seq <- NA
    }
    if(length(clones)==1 & nmuts > 0){
      i.seg <- which(colnames(seqmat)=='seg')
      i.cn <- which(colnames(seqmat)=='mutated.copies')
      i.allele <- which(colnames(seqmat)=='allele')
    }
    newclone <- list('cn'=newclone.cn, 'seq'=newclone.seq)
    clones <- c(clones, list(newclone))
  }
  if(norm.contam==TRUE){
    sampled <- c(1, sample(2:length(clones),K-1, replace=FALSE))
  }else{
    sampled <- sample(2:length(clones),K, replace=FALSE)
  }
  clones.final <- lapply(1:length(sampled),function(i){
    cndf <- as.data.frame(clones[[sampled[i]]]$cn)
    if(length(which(!is.na(unlist(clones[[sampled[i]]]$seq))))){
      seqdf <- as.data.frame(clones[[sampled[i]]]$seq)
      for(j in which(colnames(seqdf)!='allele')){
        seqdf[,j] <- as.numeric(as.character(seqdf[,j]))
        seqdf <- na.omit(seqdf[with(seqdf, order(seg, start)),])
        rownames(seqdf) <- NULL
        seqdf
      }
      seqdf$normal.copies <- na.omit(unlist(lapply(1:nrow(cndf),function(j){
        total.cn <- cndf$A[j] + cndf$B[j]
        if(length(which(seqdf$seg==j))>0){
          normal.cns <- total.cn - seqdf$mutated.copies[which(seqdf$seg==j)]
        }else{
          normal.cns <- NA
        }
        normal.cns
      })))
    }else{
      seqdf <- NA
    }
    if(nu>0){
      output <- list('cn'=cndf,'seq'=seqdf)
    }else{
      output <- list('cn'=cndf)
    }
    output
  })
  list('clones'=clones.final, 'psi'=psi)
}

snpDataGen <- function(tumor, snps.cgh=600000, sigma0.lrr=.01, sigma0.baf=.01, density.sigma=.1){
  psi <- tumor$psi
  cn.clones <- lapply(1:length(tumor$clones), function(i){tumor$clones[[i]]$cn})
  chr <- cn.clones[[1]]$chr
  eta <- data.frame('A'=Reduce('+', lapply(1:length(which(psi>0)),function(j){psi[j]*cn.clones[[j]]$A})),
                    'B'=Reduce('+', lapply(1:length(which(psi>0)),function(j){psi[j]*cn.clones[[j]]$B})))
  eta.total <- eta$A + eta$B
  mu.baf <- eta$B/eta.total
  mu.lrr <- log10(eta.total/2)
  snp.dens <- rbeta2(1:nrow(eta), .5, density.sigma)
  markers <- as.vector(rmultinom(1, size=snps.cgh, prob=snp.dens))
  sigmas.lrr <- sigma0.lrr/(markers)^.5
  sigmas.baf <- sigma0.baf/(markers)^.5
  lrr <- rnorm(length(markers), mean=mu.lrr, sd=sigmas.lrr)
  baf <- rnorm(length(markers), mean=mu.baf, sd=sigmas.baf)
  baf[baf<0] <- -baf[baf<0]
  baf[baf>1] <- 1/baf[baf>1]
  total.backcomp <- 2*(10^(lrr))
  X <- total.backcomp*(baf)
  Y <- total.backcomp*(1 - baf)
  segdf <- data.frame('chr'=chr, seg=1:nrow(eta), 'LRR'=lrr, 'BAF'=baf,'X'=X, 'Y'=Y, 'markers'=markers)
  segdf
}

#snp.info is an option for an already made snp df (like from snpDataGen) to generate
#read counts from those in addition to generating new ones.
seqDataGen <- function(tumor, snps.seq=1000000, density.sigma, mu, sigma.reads){
  psi <- tumor$psi
  cn.clones <- lapply(1:length(tumor$clones), function(i){tumor$clones[[i]]$cn})
  chrs <- cn.clones[[1]]$chr
  psi <- tumor$psi
  eta <- data.frame('A'=Reduce('+', lapply(1:length(which(psi>0)),function(j){psi[j]*cn.clones[[j]]$A})),
                    'B'=Reduce('+', lapply(1:length(which(psi>0)),function(j){psi[j]*cn.clones[[j]]$B})))
  eta.total <- eta$A + eta$B
  snp.dens <- rbeta2(1:nrow(eta), .5, density.sigma)
  markers <- as.vector(rmultinom(1, size=snps.seq, prob=snp.dens))
  starts <- cn.clones[[1]]$start
  ends <- cn.clones[[1]]$ends
  snpdf <- matrix(NA, nrow=sum(markers), ncol=7)
  colnames(snpdf) <- c('chr', 'seg', 'mut.id', 'refCounts', 'varCounts','VAF', 'totalCounts')
  for(j in 1:length(markers)){
    z <- sample(1:2, 2, replace=FALSE)
    pR <- eta[j,z[1]]/(eta[j,z[1]]+eta[j,z[2]])
    totalCounts <- round(rnorm(markers[j], 2*mu, sigma.reads))
    totalCounts[totalCounts<0] <- -totalCounts[totalCounts<0]
    refCounts <- rbinom(markers[j], totalCounts, pR)
    varCounts <- totalCounts - refCounts
    chr <- rep(chrs[j], markers[j])
    VAFs <- varCounts/(totalCounts)
    seg <- rep(j,markers[j])
    for(k in 1:length(refCounts)){
      snpdf[sum(markers[1:j]) - markers[j] +k,] <- c(chr[k], seg[k], NA, refCounts[k], varCounts[k], VAFs[k], totalCounts[k])
    }
  }
  snpdf <- as.data.frame(snpdf)
  snpdf$status <- rep('germline', nrow(snpdf))
  seq.clones <- lapply(1:length(tumor$clones), function(i){tumor$clones[[i]]$seq})
  cn.clones <- lapply(1:length(tumor$clones), function(i){tumor$clones[[i]]$cn})
  mutids <- unique(unlist(lapply(1:length(which(psi>0)),function(j){seq.clones[[j]]$mut.id})))
  mutdf <- matrix(NA, nrow=length(mutids), ncol=7)
  colnames(mutdf) <- c('chr', 'seg', 'mut.id', 'refCounts', 'varCounts','VAF', 'totalCounts')
  if(length(mutids)>0){
    for(j in 1:length(mutids)){
      indices <- which(sapply(1:length(which(psi>0)),function(k){mutids[j] %in% seq.clones[[k]]$mut.id}))
      indices.unmutated <- which(sapply(1:length(which(psi>0)),function(k){!mutids[j] %in% seq.clones[[k]]$mut.id}))
      coefs <- psi[indices]
      mutated <- sapply(indices,function(k){seq.clones[[k]]$mutated.copies[which(seq.clones[[k]]$mut.id==mutids[j])]})
      normal <- sapply(indices,function(k){seq.clones[[k]]$normal.copies[which(seq.clones[[k]]$mut.id==mutids[j])]})
      seg <- seq.clones[[indices[1]]]$seg[which(seq.clones[[indices[1]]]$mut.id==mutids[j])]
      if(length(indices.unmutated)>0){
        coefs.unmutated <- psi[indices.unmutated]
        eta.unmutated <- sum(psi[indices.unmutated]*sapply(indices.unmutated,function(i){
          cn.clones[[i]]$A[seg]+cn.clones[[i]]$B[seg]}))
      }else{
        eta.unmutated <- 0
      }
      eta <- c(sum(coefs*mutated), sum(coefs*normal) + eta.unmutated)
      pR <- eta[2]/(eta[2]+eta[1])
      totalCount <- round(rnorm(1, 2*mu, sigma.reads))
      totalCount[totalCount<0] <- -totalCount[totalCount<0]
      refCount <- rbinom(1, totalCount, pR)
      varCount <- totalCount - refCount
      VAF <- varCount/(varCount+refCount)
      chr <- seq.clones[[indices[1]]]$chr[which(seq.clones[[indices[1]]]$mut.id==mutids[j])]
      mutdf[j,] <- c(chr, seg, mutids[j], refCount, varCount, VAF, totalCount)
    }
  }
  mutdf <- as.data.frame(mutdf)
  mutdf$status <- rep('somatic', nrow(mutdf))
  vardf <- rbind(snpdf, mutdf)
  vardf <- vardf[with(vardf, order(seg)),]
  vardf
}

dataGen <- function(tumor, snps.seq, snps.cgh, mu, sigma.reads, sigma0.lrr, sigma0.baf, density.sigma){
  if(snps.cgh>0){
    cndat <- snpDataGen(tumor, snps.cgh, sigma0.lrr, sigma0.baf, density.sigma)
  }else{
    cndat <- NA
  }
  if(snps.seq){
    seqdat <- seqDataGen(tumor, snps.seq, density.sigma, mu, sigma.reads)
  }else{
    seqdat <- NA
  }
  list('cn.data'=cndat, 'seq.data'=seqdat)
}

#A plot function to visualize data and verify that the simulation is working.
plot.data <- function(tumor, data, snp=TRUE, somatic=TRUE, germline=FALSE){
  snpdata <- data$cn.data
  seqdata <- data$seq.data
  psi <- tumor$psi
  cn.clones <- lapply(1:length(tumor$clones), function(i){tumor$clones[[i]]$cn})
  psi <- tumor$psi
  markers <- snpdata$markers
  starts <- sapply(1:length(markers),function(j){sum(markers[1:j])-markers[j]})
  ends <- sapply(1:length(markers),function(j){sum(markers[1:j])})
  eta <- data.frame('A'=Reduce('+', lapply(1:length(psi),function(j){psi[j]*cn.clones[[j]]$A})),
                    'B'=Reduce('+', lapply(1:length(psi),function(j){psi[j]*cn.clones[[j]]$B})))
  lesser <- sapply(1:nrow(eta), function(j){min(eta$A[j], eta$B[j])})
  lesser.data <- sapply(1:nrow(eta), function(j){min(snpdata$X[j], snpdata$Y[j])})
  greater <- sapply(1:nrow(eta), function(j){max(eta$A[j], eta$B[j])})
  greater.data <- sapply(1:nrow(eta), function(j){max(snpdata$X[j], snpdata$Y[j])})
  errors <- c(lesser.data - lesser, greater.data - greater)
  par(mfrow=c(3,1))
  plot(0, type='n', xlim=c(0,sum(markers)), ylim=c(0,5), main='Number of Lesser Allele Copies',
       xlab='Marker Index', ylab='Expected Copy Number')
  sapply(1:nrow(eta),function(j){
    segments(x0=starts[j],x1=ends[j],y0=lesser[j],y1=lesser[j], col='black')
    segments(x0=starts[j],x1=ends[j],y0=lesser.data[j],y1=lesser.data[j], col='red')
  })
  plot(0, type='n', xlim=c(0,sum(markers)), ylim=c(0,5), main='Number of Greater Allele Copies',
       xlab='Marker Index', ylab='Expected Copy Number')
  sapply(1:nrow(eta),function(j){
    segments(x0=starts[j],x1=ends[j],y0=greater[j],y1=greater[j], col='black')
    segments(x0=starts[j],x1=ends[j],y0=greater.data[j],y1=greater.data[j], col='blue')
  })
  plot(density(errors), xlab='Error Index', main='Density of Segment Errors')
  mutations <- seqdata[seqdata$status=='somatic',]
}
mix <- function(datList, psi, datapath, simpath, index, version, alter=FALSE){
  breaks <- lapply(1:22, function(j){
    brks <- sort(unique(unlist(sapply(1:length(datList), function(k){c(datList[[k]][datList[[k]]$chrom==j,]$loc.start, 
      datList[[k]][datList[[k]]$chrom==j,]$loc.end)}))))
    starts <- brks[1:(length(brks)-1)]
    ends <- c(starts[2:length(starts)]-1, brks[length(brks)])
    data.frame('start'=starts, 'end'=ends, 'chr'=rep(j, length(starts)))
  })
  coords <- Reduce(rbind, breaks)
  columns <- sample(1:2,length(psi),replace=TRUE)
  shifts <- sapply(1:length(datList),function(j){
    lrrs <- unlist(sapply(1:nrow(datList[[j]]),function(k){rep(datList[[j]]$seg.median[k],datList[[j]]$num.mark[k])}))
    dens <- density(na.omit(2*10^lrrs))
    dens$x[which.max(dens$y)] - 2
  })
  d2 <- lapply(1:length(datList),function(j){
    temp <- datList[[j]]
    res <- t(sapply(1:nrow(coords),function(k){
      chrdat <- temp[as.character(temp$chrom)==as.character(coords$chr[k]),]
      seg <- chrdat[chrdat$loc.end>=coords$start[k] & chrdat$loc.start<=coords$end[k],]
      #Afreq <- seg$AA*2 + seg$AB
      #Bfreq <- seg$BB*2 + seg$AB
      #baf <- Bfreq/(Afreq + Bfreq)
      baf <- seg$AvgBAF
      tcn <- 2*10^(seg$seg.median)
      tcn <- tcn - shifts[j]
      if(length(which(tcn<0))>0){
        tcn[which(tcn<0)] <- -tcn[which(tcn<0)]
      }
      X <- tcn*(1-baf)
      Y <- tcn*baf
      if(nrow(seg)>0){
        output <- c(X, Y)
      }else{
        output <- rep(NA, 2)
      }
      output
    }))
    colnames(res) <- c('X', 'Y')
    res
  })
  data.mixed <- as.data.frame(Reduce('+', lapply(1:length(d2),function(k){psi[k]*d2[[k]]})))
  notna <- which(sapply(1:nrow(data.mixed),function(x){length(which(is.na(data.mixed[x,])))==0}))
  for(j in 1:length(d2)){
    d2[[j]] <- d2[[j]][notna,]
  }
  data.mixed <- data.mixed[notna,]
  starts <- coords$start[notna]
  ends <- coords$end[notna]
  chrs <- coords$chr[notna]
  data.mixed$chr <- chrs
  data.mixed$seg <- 1:length(notna)
  data.mixed$LRR <- log10((data.mixed$X + data.mixed$Y)/2)
  data.mixed$BAF <- data.mixed$X/(data.mixed$X + data.mixed$Y)
  data.mixed$markers <- unlist(sapply(1:22,function(i){
    pos.chr <- pos[pos$Chr==i,]
    data.chr <- data.mixed[data.mixed$chr==i,]
    sapply(1:nrow(data.chr),function(j){
      length(which(pos.chr >= starts[which(data.mixed$chr==i)[j]] & pos.chr <= ends[which(data.mixed$chr==i)[j]]))
    })
  }))
  data.mixed$start <- starts
  data.mixed$end <- ends
  data.mixed <- data.mixed[,c(3,4,5,6,1,2,7,8,9)]
  clones <- lapply(1:length(d2),function(k){
    #tcn <- d2[[k]][,1] + d2[[k]][,2]
    #BAF <- d2[[k]][,1]/tcn
    #BAF[which(tcn==0)] <- 0
    #df <- as.data.frame(round(t(sapply(1:length(tcn),function(l){tcn[l]*c(BAF[l], 1-BAF[l])}))))
    df <- data.frame('A'=rep(1,nrow(d2[[1]])),'B'=rep(1,nrow(d2[[1]])))
    colnames(df) <- c('A','B')
    df$chr <- data.mixed$chr
    df$start <- starts
    df$end <- ends
    df$seg <- 1:length(starts)
    df$parent.index <- rep(0, length(starts))
    df$markers <- data.mixed$markers
    df <- df[,c(3,4,5,1,2,6,7,8)]
    list('cn'=df, 'seq'=NULL)
  })
  for(k in 2:nrow(clones[[1]]$cn)){
    bin1 <- data.mixed$X[k]==data.mixed$X[k-1] & data.mixed$Y[k]==data.mixed$Y[k-1]
    bin2 <- clones[[1]]$cn$chr[k]==clones[[1]]$cn$chr[k-1]
    if(is.na(bin1)){
      bin1 <- FALSE
    }
    if(bin1==TRUE & bin2==TRUE){
      for(l in 1:length(clones)){
        clones[[l]]$cn$start[k] <- clones[[l]]$cn$start[k-1]
        clones[[l]]$cn$markers[k] <- clones[[l]]$cn$markers[k] + clones[[l]]$cn$markers[k-1]
        clones[[l]]$cn[k-1,] <- rep(NA, ncol(clones[[1]]$cn))
      }
      data.mixed$start[k] <- data.mixed$start[k-1]
      data.mixed$markers[k] <- data.mixed$markers[k-1]
      data.mixed[k-1,] <- rep(NA, ncol(data.mixed))
    }
  }
  notna <- unique(unlist(lapply(1:length(clones),function(j){which(!is.na(clones[[j]]$cn$seg))})))
  for(k in 1:length(clones)){
    clones[[k]]$cn <- clones[[k]]$cn[notna,]
    clones[[k]]$cn$seg <- rownames(clones[[k]]$cn) <- 1:nrow(clones[[k]]$cn)
    clones[[k]]$cn
  }
  data.mixed <- data.mixed[notna,]
  dat <- list('cn.data'=data.mixed, 'seq.data'=NULL)
  dat$cn.data <- na.omit(dat$cn.data)
  dat$cn.data$seg <- rownames(dat$cn.data) <- 1:nrow(dat$cn.data)
  if(alter){
    alt.pools <- lapply(1:length(psi),function(j){
      x <- dat$cn.data$X
      y <- dat$cn.data$Y
      markers <- dat$cn.data$markers
      which(markers>1000 & x-shifts[j]/2>.95 & x-shifts[j]/2<1.05 & y-shifts[j]/2>.95 & y-shifts[j]/2<1.05)
    })
    altN <- sample(1:3,length(psi),prob=c(.6,.3,.1),replace=TRUE)
    altered <- lapply(1:length(psi),function(j){sample(alt.pools[[j]],altN[j])})
    change <- lapply(1:length(psi),function(j){sample(c(-1,1),altN[j],replace=TRUE)})
    for(l in 1:length(psi)){
      if(length(altered[[l]])>0){
        for(m in 1:length(altered[[l]])){
          allele <- sample(c(1,2),1)
          clones[[l]]$cn[altered[[l]][[m]],allele+3] <- 
            clones[[l]]$cn[altered[[l]][[m]],allele+3] + change[[l]][[m]]
          dat$cn.data[altered[[l]][[m]],allele+4] <- dat$cn.data[altered[[l]][[m]],allele+4] + 
            change[[l]][[m]]*psi[l]
          if(dat$cn.data[altered[[l]][[m]],allele+4] < 0){
            dat$cn.data[altered[[l]][[m]],allele+4] <- -dat$cn.data[altered[[l]][[m]],allele+4]
          }
        }
      }
    }
  }
  tumor <- list('clones'=clones, 'psi'=psi, 'altered'=altered,'change'=change)
  save(dat, file=paste(datapath, '/mixdat', version, '-', index, '.rda', sep=''))
  save(tumor, file=paste(simpath, '/mixsim', version, '-', index, '.rda', sep=''))
}
plotSetTraits <- function(setDir, seq=TRUE, cn=TRUE, merge=FALSE, kmax=5){
  simfiles <- list.files(setDir)
  metadata <- simfiles[1]
  simfiles <- simfiles[2:length(simfiles)]
  sims <- lapply(1:length(simfiles),function(j){get(load(paste(setDir, '\\', simfiles[j], sep='')))})
  if(length(sims[[1]])==1){
    sims <- lapply(1:length(sims),function(x){sims[[x]]$tumor})
  }
  ks <- sapply(1:length(sims),function(j){length(which(sims[[j]]$psi>0))})
  psis <- t(sapply(1:length(sims),function(j){
    psi <- sims[[j]]$psi
    if(length(psi)<kmax){
      psi <- c(psi, rep(0, kmax-length(psi)))
    }
    psi
  }))
  frac.altered <- sapply(1:length(sims),function(j){
    sim <- sims[[j]]
    lengths <- sim$clones[[1]]$cn$end + 1 - sim$clones[[1]]$cn$start
    sum(lengths[unique(unlist(lapply(1:length(sim$clones),function(x){
      which(sim$clones[[x]]$cn$A!=1 | sim$clones[[x]]$cn$B!=1)})))])/sum(lengths)
  })
  if(merge){
    sim <- sims[[j]]
    lengths <- sim$clones[[1]]$cn$end + 1 - sim$clones[[1]]$cn$start
    altered <- sapply(1:length(sim$clones),function(x){sim$clones[[x]]$cn$A!=1 | sim$clones[[x]]$cn$B!=1})
    altered.indices <- which(altered)
    diffs <- diff(altered.indices)
    #To be continued?
  }else{
    length.altered <- lapply(1:length(sims),function(j){
      sim <- sims[[j]]
      lengths <- sim$clones[[1]]$cn$end + 1 - sim$clones[[1]]$cn$start
      lengths[unique(unlist(lapply(1:length(sim$clones),function(x){
        which(sim$clones[[x]]$cn$A!=1 | sim$clones[[x]]$cn$B!=1)})))]
    })
    nCNVs <- sapply(1:length(length.altered),function(z){length(length.altered[[z]])}) 
  }
  length.altered <- unlist(length.altered)
  mutations <- sapply(1:length(sims),function(j){
    sim <- sims[[j]]
    if(length(sim)==1){
      sim <- sim$tumor
    }
    length(unique(unlist(lapply(1:length(sim$clones), function(x){sim$clones[[x]]$seq$mut.id}))))
  })
  if(cn==TRUE & seq==FALSE){
    #par(mfrow=c(1,2))
    #hist(ks, breaks=c(.5,1.5,2.5,3.5,4.5,5.5), xlab='Number of Clones', main='')
    plot(density(unlist(psis),bw=.05), xlab='Psi', xlim=c(0,1), main='F) Density of Psi Values across Simulations')
    plot(density(frac.altered,bw=.01), xlab='Genome Altered Fraction', xlim=c(0,1), main='G) Genome AF Density among Simulations')
    #plot(density(length.altered,bw=.01), xlab='CNV Lengths', xlim=c(0,max(length.altered)), main='')
    #hist(nCNVs, xlab='CNVs per Sample', main='')
  }else if(cn==FALSE & seq==TRUE){
    #par(mfrow=c(1,2))
    #hist(ks, breaks=c(.5,1.5,2.5,3.5,4.5,5.5), xlab='Number of Clones', main='')
    plot(density(unlist(psis),bw=.05), xlab='Psi', xlim=c(0,1), main='A) Density of Psi Values across Simulations')
    plot(density(mutations,bw=5), xlab='Mutations per Sample',xlim=c(0,max(mutations)), main='B) Density of Mutations per Simulation')
  }else if(cn & seq){
    #par(mfrow=c(1,3))
    #hist(ks, breaks=c(.5,1.5,2.5,3.5,4.5,5.5), xlab='Number of Clones', main='')
    plot(density(unlist(psis),bw=.05), xlab='Psi', xlim=c(0,1), main='C) Density of Psi Values across Simulations')
    plot(density(frac.altered,bw=.01), xlab='Genome Altered Fraction', xlim=c(0,1), main='D) Genome AF Density among Simulations')
    #plot(density(length.altered,bw=.01), xlab='CNV Lengths', xlim=c(0,max(length.altered)), main='')
    #hist(nCNVs, xlab='CNVs per Sample', main='')
    plot(density(mutations,bw=5), xlab='Mutations per Sample',xlim=c(0,max(mutations)), main='E) Density of Mutations per Simulation')
  }
}
plotSNP <- function(chr.sim, chr.dat, j, sam, homFrac=.68, mu.0=0, sigma.0=.03){
  dat <- get(load(paste('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\dat\\data-3', '\\dat-3-', j,'.rda', sep='')))
  sim <- get(load(paste('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-3', '\\sim-3-', j,'.rda', sep='')))
  segdata <- dat$dat$cn.data
  #sigma0.lrr <- metadata$data.params$sigma0.lrr
  sigma0.lrr <- .12
  sigma0.baf <- metadata$data.params$sigma0.baf
  segs <- segdata[segdata$chr==chr.sim,]
  x.mu <- unlist(sapply(1:nrow(segs),function(k){rep(segs$X[k], segs$markers[k])}))
  y.mu <- unlist(sapply(1:nrow(segs),function(k){rep(segs$Y[k], segs$markers[k])}))
  index <- 1:length(x.mu)
  hom <- sample(1:length(index), round(length(index)*homFrac), replace=FALSE)
  inverted <- sample(1:length(index), round(length(index)*.5), replace=FALSE)
  x.mu[hom] <- x.mu[hom] + y.mu[hom]
  y.mu[hom] <- 0
  x.inverted <- x.mu[inverted]
  x.mu[inverted] <- y.mu[inverted]
  y.mu[inverted] <- x.inverted
  lrr.mu <- log10((x.mu + y.mu)/2)
  baf.mu <- x.mu/(x.mu + y.mu)
  lrr <- rnorm(length(lrr.mu), lrr.mu, sigma0.lrr)
  baf <- rnorm(length(baf.mu), baf.mu, sigma0.baf)
  baf[which(baf<0)] <- -baf[which(baf<0)]
  baf[which(baf>1)] <- 1/baf[which(baf>1)]
  path <- paste('E:\\SNP-HapMap\\MergeByChrom\\Full Data Table_Chr_', chr.dat, '.txt', sep='')
  cndata <- read.table(path, header=T)
  columns <- unique(sapply(4:ncol(cndata),function(x){strsplit(colnames(cndata)[x],split='[.]')[[1]][1]}))
  samn <- which(columns==sam)
  lrr.dat <- cndata[,samn*3+2]
  baf.dat <- cndata[,samn*3+3]
  par(mfrow=c(2,2))
  #maintab <- paste('Chr ', chr, sep='')
  plot(1:length(lrr), lrr, xlab='SNP Index', ylab='LRR', ylim=c(-1.5,1), pch='.', col='gray',
       main='A) SNP Array Data (simulated data)')
  plot(1:length(lrr.dat), lrr.dat, xlab='SNP Index', ylab='LRR', ylim=c(-1.5,1), pch='.', col='gray',
       main='B) SNP Array Data (real data)')
  abline(h=c(log10((1:4)/2)), col=c('gray','black','gray','gray'))
  abline(h=c(log10((1:4)/2)), col=c('gray','black','gray','gray'))
  plot(1:length(baf), baf, main='', xlab='SNP Index', ylab='BAF', ylim=c(0,1), pch='.', col='gray')
  abline(h=.5, col='black')
  plot(1:length(baf.dat), baf.dat, main='', xlab='SNP Index', ylab='BAF', ylim=c(0,1), pch='.', col='gray')
  abline(h=.5, col='black')
  #plot(density(baf), xlim=c(-1.5,1.5), ylim=c(0,5), main=paste('SD=',round(sd(baf),digits=3), sep=''))
}

plotSEQ <- function(seqdata, seqsim){
  som.dat <- seqdata[seqdata$somatic_status=='SOMATIC',]
  vaf.dat <- na.omit(som.dat$alt_read_count/(som.dat$alt_read_count + som.dat$ref_read_count))
  tc.dat <- som.dat$alt_read_count + som.dat$ref_read_count
  seq.sim <- seqsim$dat$seq.data
  som.sim <- seq.sim[seq.sim$status=='somatic',]
  tc.sim <- som.sim$totalCounts
  vaf.sim <- som.sim$VAF
  par(mfrow=c(2,2))
  plot(tc.sim[1:225], ylim=c(0,250), xlab='Index', ylab='Total Read Counts', main='A) Somatic Variants (Simulation)', pch=16)
  plot(tc.dat[1:225], ylim=c(0, 250), xlab='Index', ylab='Total Read Counts', main='B) Somatic Variants (Data)', pch=16)
  plot(vaf.sim[1:225], ylim=c(0, 1), xlab='Index', ylab='VAF', main='', pch=16)
  plot(vaf.dat[1:225], ylim=c(0, 1), xlab='Index', ylab='VAF', main='', pch=16)
}

doSimGen <- function() {
  library(gtools)
source('C:\\Users\\Mark\\OneDrive - The Ohio State University\\Rcode\\functions.R')
chrtab <- read.table('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sciClone\\chlens.txt')
ord <- sort(as.numeric(as.character(chrtab[,1])), index.return=TRUE)$ix
chrtab <- chrtab[ord,]
chlens <- as.numeric(chrtab[,2])
rm(eval)
#source("https://bioconductor.org/biocLite.R")
#biocLite("DNAcopy")
#library(DNAcopy)
###Generate a sample simulation population:
#set: 1 (pcnv=0, nu=30, sigma.lrr=.15, sigma.baf=.03, sigma.reads=25); k = 1:5
#set: 2 (pcnv=1, nu=30, sigma.lrr=.15, sigma.baf=.03, sigma.reads=25); k = 1:5
#set: 3 (pcnv=1, nu=0, sigma.lrr=.15, sigma.baf=.03, sigma.reads=25); k = 1:5
gen.sims <- FALSE

#path1 <- 'C:\\Users\\Mark\\OneDrive - The Ohio State University\\Integrated\\parOpt\\sim'
#path2 <- 'C:\\Users\\Mark\\OneDrive - The Ohio State University\\Integrated\\parOpt\\data'
path1 <- 'C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim2'
path2 <- 'C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\dat2'
if(gen.sims){
  pars.default <- list('rounds'=400, 'nu'=0, 'pcnv'=1, 'norm.contam'=FALSE)
  dataPars <- list('snps.seq'=1000000,'snps.cgh'=600000,'mu'=70,'sigma.reads'=25,'sigma0.lrr'=.15, 
                   'sigma0.baf'=.03,'density.sigma'=.1)
  threshold <- .04
  dirs <- list.dirs(path1)
  dirs <- dirs[-1]
  if(length(dirs)>0){
    setIDs <- as.numeric(sapply(1:length(dirs),function(j){strsplit(dirs[j], split='sims-')[[1]][2]}))
    setID <- max(setIDs) + 1
  }else{
    setID <- 1
  }
  dir.create(paste(path1, '\\sims-', setID, sep=''))
  dir.create(paste(path2, '\\data-', setID, sep=''))
  simpath <- paste(path1, '\\sims-', setID, sep='')
  datapath <- paste(path2, '\\data-', setID, sep='')
  #psis <- genSimplex(20, 5, reps=1)
  psis <- get(load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\Integrated\\psis.100.rda'))
  psis <- psis[which(sapply(1:nrow(psis),function(j){length(which(psis[j,]==0 | psis[j,]>threshold))==ncol(psis)})),]
  ks <- sapply(1:nrow(psis),function(j){length(which(psis[j,]>0))})
  seed <- 123
  set.seed(seed)
  sampled <- unlist(lapply(1:5, function(x){sample(which(ks==x), 60, replace=TRUE)}))
  psis <- psis[sampled,]
  n <- nrow(psis)
  for(j in 1:nrow(psis)){
    go <- FALSE
    while(go==FALSE){
      tumor <- try(tumorGen(psis[j,], pars.default$rounds, pars.default$nu, pars.default$pcnv, pars.default$norm.contam),silent=TRUE)
      if(class(tumor)!='try-error'){
        data <- dataGen(tumor, dataPars$snps.seq, dataPars$snps.cgh, dataPars$mu, dataPars$sigma.reads, 
                        dataPars$sigma0.lrr, dataPars$sigma0.baf, dataPars$density.sigma)
      }
      if(class(tumor)!='try-error' & class(data)!='try-error'){
        go <- TRUE
      }
    }
    filename.tumor <- paste(simpath, '\\', 'sim-', setID, '-', j, '.rda', sep='')
    filename.data <- paste(datapath, '\\', 'dat-', setID, '-', j, '.rda', sep='')
    sim <- list('tumor'=tumor)
    save(sim, file=filename.tumor)
    dat <- list('dat'=data)
    save(dat, file=filename.data)
  }
  if(pars.default$nu > 0 & pars.default$pcnv > 0){
    alt.types <- 'both'
  }else if(pars.default$nu > 0 & pars.default$pcnv == 0){
    alt.types <- 'mut'
  }else if(pars.default$nu == 0 & pars.default$pcnv > 0){
    alt.types <- 'cnv'
  }else{
    alt.types <- 'none'
  }
  metadata <- list('tumor.params'=pars.default, 'data.params'=dataPars, 'alt.types'=alt.types, 'date'=Sys.Date(),
                  'setID'=setID, 'seed'=seed)
  save(metadata, file=paste(simpath, '\\metadata-', setID, '.rda', sep=''))
}
gen.mixtures <- FALSE
if(gen.mixtures){
  version <- 1
  threshold <- .04
  if(version==1){
    version <- ''
  }
  seed <- 12345
  set.seed(seed)
  path <- 'C:/Users/Mark/OneDrive - The Ohio State University/clonetools'
  rewind.hapmap <- get(load(paste(path, '\\SNP_array_data\\rewind.hapmap.rda', sep='')))
  rewind.cll <- get(load(paste(path, '\\', 'SNP_array_data\\rewind_both.rda', sep='')))
  rewind.cll <- rewind.cll[,c(1:6,11)]
  rewind <- rbind(rewind.cll, rewind.hapmap)
  datapath <- paste('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/dat/data-mix2',version,sep='')
  simpath <- paste('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/sim/sims-mix2',version,sep='')
  psipool <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/clonefinder/psis.100.rda'))
  psipool <- psipool[which(sapply(1:nrow(psipool),function(j){length(which(psipool[j,]==0 | psipool[j,]>threshold))==ncol(psipool)})),]
  ks <- sapply(1:nrow(psipool),function(i){length(which(psipool[i,]>0))})
  N <- 300 #Let's to 60 apiece for k=1, 2, 3, 4, and 5
  indices <- c(sample(which(ks==1),60,replace=TRUE),sample(which(ks==2),60,replace=TRUE),sample(which(ks==3),60,replace=FALSE),
               sample(which(ks==4),60,replace=FALSE),sample(which(ks==5),60,replace=FALSE))
  psis <- lapply(1:N,function(j){psipool[indices[j],][which(psipool[indices[j],]>0)]})
  #IDpool <- c('CL001','CL060','CL082','CL084','CL102','CL185','NA19193')
  IDpool <- c('NA07019', 'NA12234', 'NA12249', 'NA12753', 'NA12761', 'NA18545', 'NA18975', 'NA18999', 'NA18517')
  IDs <- lapply(1:N, function(j){sample(IDpool, length(which(psis[[j]]>0)), replace=FALSE)})
  metadata <- list('IDs'=IDs, 'psis'=psis, 'seed'=seed)
  save(metadata, file=paste(simpath, '\\metadata.rda', sep=''))
  pos <- get(load('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/pos.rda'))
  sapply(1:length(psis),function(i){mix(lapply(1:length(IDs[[i]]),function(j){
    rewind[rewind$SamID==IDs[[i]][j],]}),psis[[i]], datapath, simpath, i, version, alter=TRUE)})
}

check <- FALSE
if(check){
  load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-mix\\metadata-mix.rda')
  index <- 196
  sams <- metadata$IDs[[index]]
  first.two1 <- paste(strsplit(sams[1],split='')[[1]][1:2],collapse='')
  first.two2 <- paste(strsplit(sams[2],split='')[[1]][1:2],collapse='')
  folder1 <- folder2 <- 'data-CLL'
  if(first.two1=='NA'){
    folder1 <- 'data-hapmap'
  }
  if(first.two2=='NA'){
    folder2 <- 'data-hapmap'
  }
  dat1 <- get(load(paste('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\dat\\',
                         folder1,'\\dat-',sams[1],'.rda', sep='')))
  dat2 <- get(load(paste('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\dat\\',
                         folder2,'\\dat-',sams[2],'.rda', sep='')))
  mix.dat <- get(load(paste('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\dat\\data-mix\\mixdat-', index, '.rda', sep='')))
  dat1 <- dat1$dat$cn.data
  dat2 <- dat2$dat$cn.data
  mix.dat <- mix.dat$cn.data
  psi <- metadata$samplePsis[index,]
  cnplotfun(dat1, 2, 'dat1')
  cnplotfun(dat2, 2, 'dat1')
  cnplotfun(mix.dat, 2, 'mixture')
  mixplotfun(mix.dat, dat1, dat2, 2)
}

png('C:/Users/Mark/OneDrive - The Ohio State University/clonetools/cf_images/fig_2.png',height=500,width=900,res=100)
layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5,6,6,6,7,7,7), 3, 6, byrow = TRUE))
plotSetTraits('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-1', cn=FALSE)
plotSetTraits('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-2')
plotSetTraits('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-3', seq=FALSE)
dev.off()

#png('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\images\\sims-1.png', width=10, height=3.5, units="in", res=300)
#plotSetTraits('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-1', cn=FALSE)
#dev.off()
#png('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\images\\sims-2.png', width=15, height=3.5, units="in", res=300)
#plotSetTraits('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-2')
#dev.off()
#png('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\images\\sims-3.png', width=10, height=3.5, units="in", res=300)
#plotSetTraits('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-3', seq=FALSE)
#dev.off()
#png('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\images\\sims-mix.png', width=10, height=3.5, units="in", res=300)
#plotSetTraits('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-mix', seq=FALSE, merge=FALSE)
#dev.off()

###
imgpath <- 'C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\images'
metadata <- get(load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\sim\\sims-3\\metadata-3.rda'))
plotSNP(chr.sim=2, j=68, homFrac=.68, mu.0=0, sigma.0=.03, chr.dat=2, sam='NA18855')
plotSEQ(seqdata=get(load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\deepcnv\\GS1115.rda')),
    seqsim=get(load('C:\\Users\\Mark\\OneDrive - The Ohio State University\\clonetools\\dat\\data-1\\dat-1-10.rda')))

}
