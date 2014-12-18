library(MASS)
library(LearnBayes)
setwd("C:\\Users\\Mark\\Documents\\Lab\\New Data Sets")
dcll<-read.table("params_CLL.txt", header=T)
cll<-dcll[dcll$num.mark>=200,]
cll$LRR<-cll[[6]]
cll$BAF<-cll[[11]]

LOH<-cll[cll$Copies=="LOH",]
plot(LOH$LRR, LOH$BAF)
LOH<-LOH[LOH$BAF<=0.25,]
LOH.BAF<-fitdistr(LOH$BAF, "beta", list(shape1=0.2, shape2=0.1))
LOH.LRR<-fitdistr(LOH$LRR, "normal")
LOH.LRR1<-LOH.LRR[[1]]
LOH.BAF1<-LOH.BAF[[1]]

N01<-cll[cll$Copies=="N01",]
plot(N01$LRR, N01$BAF)
N01<-N01[N01$BAF<=0.25,]
N01.BAF<-fitdistr(N01$BAF, "beta", list(shape1=0.2, shape2=0.1))
N01.LRR<-fitdistr(N01$LRR, "normal")
N01.LRR1<-N01.LRR[[1]]
N01.BAF1<-N01.BAF[[1]]

N02<-cll[cll$Copies=="N02",]
plot(N02$LRR, N02$BAF)
N02<-N02[N02$BAF>=0,]
N02.BAF<-fitdistr(N02$BAF, "beta", list(shape1=0.2, shape2=0.1))
N02.LRR<-fitdistr(N02$LRR, "normal")
N02.LRR1<-N02.LRR[[1]]
N02.BAF1<-N02.BAF[[1]]

N03<-cll[cll$Copies=="N03",]
plot(N03$LRR, N03$BAF)
N03<-N03[N03$BAF>=0,]
N03.BAF<-fitdistr(N03$BAF, "beta", list(shape1=0.2, shape2=0.1))
N03.LRR<-fitdistr(N03$LRR, "normal")
N03.LRR1<-N03.LRR[[1]]
N03.BAF1<-N03.BAF[[1]]

N04<-cll[cll$Copies=="N04",]
plot(N04$LRR, N04$BAF)
N04<-N04[N04$BAF>=0,]
N04.BAF<-fitdistr(N04$BAF, "beta", list(shape1=0.2, shape2=0.1))
N04.LRR<-fitdistr(N04$LRR, "normal")
N04.LRR1<-N04.LRR[[1]]
N04.BAF1<-N04.BAF[[1]]

df<-data.frame(rbind(LOH.LRR1, N01.LRR1, N02.LRR1, N03.LRR1, N04.LRR1), rbind(LOH.BAF1, N01.BAF1, N02.BAF1, N03.BAF1, N04.BAF1))

#Sample function for generating pdfs for the LOH compartment coefficient:
FLOH<-function(x, a_y, b_y){
m<-1/2.718^((df[1, 1]-x)^2)
s<-df[1, 2]
a<-1/(2.718^((df[1, 3]-a_y)^2)*b_y)
b<-1/(2.718^((df[1, 4]-b_y)^2)*b_y)
c(m, s, a, b)
}
dLOHc1<-mapply(FLOH, cll$LRR, cll$a, cll$b)
paramsLOHc1<-data.frame(cbind(dLOHc1[1,], dLOHc1[2,], dLOHc1[3,], dLOHc1[4,]))
paramsLOHc1$segment<-cll$segment
d<-paramsLOHc1
func<-function(x, q){
dnorm(q, d[x, 1], d[x, 2])*dbeta(q, d[x, 3], d[x, 4])
}

#func2 plots the pdf at an interval of .0025 along the X-axis 
X=seq(0, 1, by=.0025)
func2<-function(x){
list=rep(x, 401)
Y=mapply(func, list, X)
plot(X, Y, main=paste("Segment ", d$segment[[x]], sep=""))
}

#Example: func2(19) shows a pdf of the LOH coefficient at 11.13.4 (segments are out of order in the data frame because they were ordered by call at 
#some point.There peak is around .8.
func2(19)
