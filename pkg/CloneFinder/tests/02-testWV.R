library(CloneFinder)

new("WeightVector")
new("WeightVector", 1:3)
try( new("WeightVector", 0) )
try( new("WeightVector", -1:3) )

WeightVector(1:4)
try( WeightVector(0) )
try( WeightVector(-1:3) )
try( WeightVector(LETTERS[1:3]) )
try( WeightVector(c(1, 2, 'a')) )

wv <- WeightVector(4:1)
as(wv, "numeric")

try( as.numeric(wv) ) # doesn't play well in the methods sandbox

canCoerce(wv, "numeric")
canCoerce(wv, "double") # more methods weirdness. "double" doesn't have a class def
