\name{Clone-class}
\Rdversion{1.1}
\docType{class}
\alias{Clone-class}
\alias{Clone}

\title{Class \code{"Clone"}}
\description{
The \code{Clone} class is used to simulate segmented copy number data
from a sample that contains a mixture of several subclones.
}
\usage{
Clone(nSegments, weights = rep(1/5, 5), segnames = NULL)
}
\arguments{
  \item{nSegments}{The number of segments.}
  \item{weights}{A numeric vector, the prevalence of each pure
    'compartment'.  Each entry is requirted to be non-negative. If
    necessary, the weights will be adjusted to ensure that they sum to
    1.  The length of this vector implicitly defines the number of
    compartments.}
  \item{segnames}{An optional character vector with names for the segments.}
}
\details{
  A "compartment" represents a collection of cells that all have the
  same pure copy number state.  For example, the default state of most
  genomic segments is to contain two different (i.e., heterozygous)
  copies.  Another possible "compartment" occurs when all cells have
  gained one copy of a chromosome.

  A "Clone" represents a homogeneous collection of cells.  In its
  implementation here, a "Clone" consists of an ordered list of
  segments, where each segment represents exactly one pure compartment.

  When you create an object ofthe "Clone" class using the \code{Clone}
  constructor, each segment is assigned to one of the pure compartments,
  with probability determined by the \code{weights} (or frequiencies). 
}
\value{
  The \code{Clone} function returns an object of the "Clone" class.
}

\section{Objects from the Class}{
Objects should be created by calls to the constructor function, \code{Clone}.
}
\section{Slots}{
  \describe{
    \item{\code{weights}:}{An bject of class \code{"numeric"}.  The
      length of the \code{weoights} vector is the number of components. }
    \item{\code{segments}:}{An integer vector.  Each entry is an integer
      between 1 and \code{length(weights)} and should be interpreted as
      an index identifying the compartment.
  }
}
\section{Methods}{
No methods are defined with class "Clone" in the signature.
}
\references{
}
\author{Kevin R. Coombes <krc@silicovore.com>,
  Mark Zucker 
}
\note{
}
\seealso{
}
\examples{
showClass("Clone")
}
\keyword{classes}
