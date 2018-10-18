\name{Simulating Clones}
\alias{generateTumorData}
\alias{plotTumorData}
\alias{tumorGen}
\alias{dataGen}
\title{ Simulating Tumor Clones }
\description{
  Simulating copy number segmentation data and sequencing mutation data
  for tumors composed of multiple clones.
}
\usage{
generateTumorData(tumor, snps.seq, snps.cgh, mu, sigma.reads,
                  sigma0.lrr, sigma0.baf, density.sigma)
plotTumorData(tumor, data, snp = TRUE, somatic = TRUE, germline = FALSE)
tumorGen(...)
dataGen(tumor, ...)
}
\arguments{
  \item{tumor}{ an object of the \code{Tumor} class. }
  \item{snps.seq}{TODO}
  \item{snps.cgh}{TODO}
  \item{mu}{TODO}
  \item{sigma.reads}{TODO}
  \item{sigma0.lrr}{TODO}
  \item{sigma0.baf}{TODO}
  \item{density.sigma}{TODO}
  \item{data}{TODO}
  \item{snp}{logical}
  \item{somatic}{logical}
  \item{germline}{logical}
  \item{...}{additional variables}
}
\details{
  TODO

  Note that \code{tumorGen} and \code{dataGen} are DEPRECATED.
}
\value{
  TODO
}
\author{Kevin R. Coombes <krc@silicovore.com>, Mark Zucker <zucker.64@buckeyemail.osu.edu>}
\examples{
# TODO
}
\keyword{array}
