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
plotTumorData(tumor, data)
tumorGen(...)
dataGen(tumor, ...)
}
\arguments{
  \item{tumor}{ an object of the \code{Tumor} class. }
  \item{snps.seq}{an integer; the total number of germline variantts and
  somatic mutations to simulate int he tumor genome.}
  \item{snps.cgh}{an integer; the number of single nucleotide
    polymorphisms (SNPs) to simulate as meaurements made to estimate
    copy number.}
  \item{mu}{an integer; the average read depth of a simulated sequencing
    study giving rise to mutations.}
  \item{sigma.reads}{a real number; the standard deviation of the number
    of simulated sequencing reads per base.}
  \item{sigma0.lrr}{a real number; the standard deviation of the
    simulated per-SNP log R ratio (LRR) for assessing copy number.}
  \item{sigma0.baf}{a real number; the standard deviation of the
    simulated B allele frequency (BAF) for assessing copy number.}
  \item{density.sigma}{a real number; the standard deviation of a beta
    distribution used to simulate the number of SNP markers per copy
    number segment.}
  \item{data}{ a list containing two data frames, \code{cn.data} and
    \code{seq.data}, as produced by \code{generateTumorData}.}
  \item{...}{additional variables}
}
\details{
  Copy number and mutation data are simulated essentially
  independently. TODO: MARK TO EXPLAIN SIMULATION MODELS.

  Note that \code{tumorGen} (an alias for \code{Tumor} that returns a
  list instead of a Tumor object) and \code{dataGen} (an alias for
  \code{generateTumorData}) are DEPRECATED.
}
\value{
  The \code{generateTumorData} function returns a list with two
  components, \code{cn.data} and \code{seq.data}. Each component is
  itself a data frame. Note that in some cases, one of these data frames
  may have zero rowes or may be returned as an \code{NA}.

  The \code{cn.data} component contains seven columns: \code{chr}, the
  chromosome number; \code{seq}, a unique segmwent identifier;
  \code{LRR}, simulate sebment-wise log ratios; \code{BAF}, simulated
  segment-wise B allele frequencies; \code{X} and \code{Y}, simulated
  intensities for two separate alleles/haplotypes per segment; and
  \code{markers}, the simulated number of SNPS per segment.

  The \code{seq.data} component contains eight columns: \code{chr} and
  \code{seq}, as above; \code{mut.id}, a unique mutation identifier;
  \code{refCounts} and \code{varCounts}, the simulated numbers of
  reference and variant counts per; \code{VAF}, the simulated variant
  allele frequency; \code{totalCounts}, the simulated total number of
  read counts; and \code{status}, a character (that should probbaly bee
  a factor) indicating whether a variant should be viewed as somatic or
  germline.

  The \code{plotTumorData} fuinction invisibly returns its \code{data}
  argument.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  Mark Zucker \email{zucker.64@buckeyemail.osu.edu}
}
\examples{
psis <- c(0.6, 0.3, 0.1) # three clones
# create tumor with copy number but no mutation data
tumor <- Tumor(psis, rounds = 400, nu = 0, pcnv = 1, norm.contam = FALSE)
# simulate the dataset
dataset <- generateTumorData(tumor, 10000, 600000, 70, 25, 0.15, 0.03, 0.1)
#plot it
plotTumorData(tumor, dataset)
}
\keyword{array}
