library(CodeDepends)
library(Rgraphviz)

home = getwd()
myenv <- new.env()

setwd("../pkg/CloneFinder/R")
source("00-generics.R", local=myenv)
source("01-cloneMaker.R", local=myenv)
source("02-prefit.R", local=myenv)
source("03-psi.R", local=myenv)
setwd(home)

setwd("TempCode")
source("functions.R", local=myenv)
source("algorithm.R", local=myenv)
source("assessing.R", local=myenv)
source("running.R", local=myenv)
source("simgen.R", local=myenv)
source("doOldTests.R", local=myenv)
setwd(home)

syms <- objects(myenv)
objs = lapply(syms, function(x) {
  f = get(x, myenv)
  if (is.function(f)) 
    f
  else NULL
})
names(objs) = syms
funs <- objs[!sapply(objs, is.null)]

mycg <- makeCallGraph(unlist(funs, recursive = FALSE), 
                      all = FALSE, recursive = TRUE,
                      names = unlist(lapply(funs, names)),
                      package="CF")
gg <- layoutGraph(mycg, layoutType="fdp")
graph.par(list(nodes=list(fontsize=40)))
renderGraph(gg)

write.table(names(mycg@edgeData), file="edgeSet.txt",
            col.names=FALSE, row.names=FALSE, quote=FALSE)


### looking for duplicate names
fudge <- read.table("../../../rforge/clonefinder/Flow/funclist.tsv", sep="\t", quote="", header=FALSE)
tab <- table(fudge$V2)
tab[tab>1]
