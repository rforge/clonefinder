library(Polychrome)

p34 <- palette36.colors(36)[3:36]
names(p34) <- colorNames(p34)
p34.deut <- colorDeficit(p34, "deut")
p34.prot <- colorDeficit(p34, "prot")
p34.trit <- colorDeficit(p34, "trit")

shift <- function(i, k=34) c(i, 1:(i-1), (1+i):k)
co <- shift(13)
pd <- computeDistances(p34.deut[co])
pp <- computeDistances(p34.prot[co])
pt <- computeDistances(p34.trit[co])

rd <- rank(pd)[order(names(pd))]
rp <- rank(pd)[order(names(pp))]
rt <- rank(pd)[order(names(pt))]
score <- 2*rd + 1.5*rp + rt

x <- p34[names(rev(sort(score)))][1:12]
y <- colorDeficit(x, "deut")
z <- colorDeficit(x, "prot")
w <- colorDeficit(x, "trit")

opar <- par(mfrow=c(2,2))
swatch(x, main="Normal")
swatch(y, main="Deuteranope")
swatch(z, main="Protanope")
swatch(w, main="Tritanope")
