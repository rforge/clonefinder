\name{Simplices}
\alias{sampleSimplex}
\alias{generateSimplex}

\title{ Simplices and Clonal Fractions }
\description{
  Utility functions for working with vectors of clonal fractions.
}
\usage{
sampleSimplex(n, d = 5)
generateSimplex(k, d, reps = 1)
}
\arguments{
  \item{d}{ an integer, the dimension of the simplex, or the number of
    clones. }
  \item{n}{ an integer, the number of vectors to sample randomly.  }
  \item{k}{ an integer, the number of equally spaced points to select along
  each side of the simplex while constructing a lattice. }
  \item{reps}{ an integer, the number of times to repeat the lattice. }
}
\details{
  When studying the clonal subpopulations of a tumor sample, we
  frequently need access to vectors that contain the fraction of cells
  belonging to each clone. These vectors are characterized by the fact
  that each entry is nonzero and they must add up to 1. The set of such
  vectors/points in d-dimensional space defines the "d-simplex". The
  functions defined here allow us to work with d-simplices, either by
  randomly sampling vectors (\code{sampleSimplex}) or by systematically
  filling the space with a regular lattice (\code{generateSimplex}).
}
\value{
  Both functions return a matrix with \code{d} columns. Each row
  contains nonzeo real numbers that sum to 1. The \code{generateSimplex}
  function ensures that (a) each row is unique and (b) the entries in
  each row appear in decreasing order. 
}
\author{Kevin R. Coombes <krc@silicovore.com>, Mark Zucker <zucker.64@buckeyemail.osu.edu>}
\examples{
sampleSimplex(5, 3)
generateSimplex(5, 3)
}
\keyword{array}
