library(CloneSeeker)

RNGversion("3.5.3")
set.seed(461283)
sampleSimplex(4, 3)
all( apply(sampleSimplex(4, 3), 1, sum) == 1) # assertion

generateSimplex(5, 3)

# These should all fail, for obvious reasons.
try( sampleSimplex(0) )
try( sampleSimplex(-1) )
try( sampleSimplex(4, 0) )
try( sampleSimplex(4, -1) )
