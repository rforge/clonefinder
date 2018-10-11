\name{WeightVector-class}
\alias{WeightVector-class}
\alias{WeightVector}
\alias{coerce,WeightVector,numeric-method}
\docType{class}
\title{Class "WeightVector"}
\description{
  A class, with validity checking, to contain vectors of non-negative
  real number whose sum equals one.
}
\usage{
WeightVector(phi)
}
\arguments{
  \item{phi}{ a numeric vector containing non-negative values.}
}
\details{
  When trying to simulate or model biological data from (sub)clonal
  populations of cells, we need vectors that keep track of the fraction
  of cells belonging to each clone. These vectors can only contain
  non-negative entries, and the entries must add up to one. (We
  thought about calling these things "ClonalFractions", but that seems
  overly specialized for a notion that is likely to prove useful in
  other contexts.) Such vectors of length \code{d} can also be viewed as
  points of a d-dimensional simplex.

  We have implemented \code{WeightVector}s as an S4 class, primarily so
  we can enforce the defining properties. We also expect this design to
  make it easier to use them as slots in other classes.
}
\section{Objects from the Class}{
  Although objects can be created using \code{new}, the preferred method
  is to use the constructor function, \code{WeightVector}.
}
\section{Slots}{
  \describe{
    \item{\code{psi}}{ a vector containing non-negative values whose sum
      equals one. }
  }
}
\section{Methods}{
  \describe{
    \item{coerce(from, to, strict = TRUE)}{Convert the
      \code{WeightVector} object into a simple numeric vector. Never
      actually used in this form, since the preferred method is to write
      \code{as(WV, "numeric")}.}
  }
}
\value{
  The constructor returns a valid object of the \code{WeightVector}
  class.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
}
\seealso{
  \code{\link{Simplices}}
}
\examples{
showClass("WeightVector")

WeightVector(1:4)
as(WeightVector(c(2, 3, 5)), "numeric")
}
\keyword{classes}

