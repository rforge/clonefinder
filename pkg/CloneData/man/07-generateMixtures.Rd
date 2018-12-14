\name{Generating data from artificial mixtures}
\alias{generateMixtures}
\title{ Generating sets of artificially mixed and altered heterogeneous data}
\description{
  Generating and saving a 'simulated' tumor data set by artificially
  mixing and altering real SNP array data that can be used in clonal
  heterogeneity analysis to assess accuracy of algorithms.
}
\usage{
generateMixtures(dataPath, mixPath, nPerK, segmentedData, ID_pool, pos)
}
\arguments{
  \item{dataPath}{ path to which simulated tumors will be saved.}
  \item{mixPath}{ path to which artificially mixed and altered SNP array
    data will be saved.}
  \item{nPerK}{a vector of integers denoting the number of tumors to
    generate for each possible number of clones, where the nth entry
    dictates how many n-clone tumors will be generated.}
  \item{segmentedData}{segmented SNP array data from which mixtures will
    be generated; must contain following columns: 'loc.start' (segment
    start locus), 'loc.end' (segment end locus), 'seg.median' (median
    Log R ratio), 'SamID' (sample ID), 'chrom' (chromosome number),
    'AvgBAF' (average B allele fraction for segment), 'num.mark' (number
    of markers per segment).}
  \item{ID_pool}{a list of sample IDs from segmentedData from which
    samples will be drawn to generate artificial mixtures.}
  \item{pos}{a data frame with two columns, \code{Chr} and
    \code{Position}, defining the chromsoomal locations of the simulated
    SNPs.}
}
\details{
  A set of artificial mixtures (with CNVs artificially added) can be
  generated from real SNP array data. The number of artificial mixtures
  to generate - and how many mixtures for each possible number of clones
  to generate - can be set with the input parameters.
}
\value{
  The \code{generateMixtures} function generates and saves two lists for
  each mixture: a 'tumor' (consisingt of artificially altered real data
  making up the 'clones' of the mixture, saved in the path 'simpath'),
  with objects: \code{psi}, a vector of clonal fractions, \code{clones},
  which is a list of tumor clones, each of which in turn consists of a
  data frame \code{cn} and a data frame \code{seq}, a list
  \code{altered} (a list of segments artificially altered), and a list
  \code{change} (the copy number change introduced to the altered
  segments); and a simulated data object (saved in the path 'datapath'),
  with objects: \code{cn.data} and \code{se .data}. Each component is
  itself a data frame. Note that in some cases, one of these data frames
  may have zero rows or may be returned as an \code{NA}.

  Each list in the \code{cn} component contains seven columns:
  \describe{
    \item{\code{chr}}{the chromosome number;}
    \item{\code{start}}{the starting locus of each genomic segment;}
    \item{\code{end}}{the ending locus of each genomic segment;}
    \item{\code{A}}{the first allelic copy number at each genomic segment;}
    \item{\code{B}}{the second allelic copy number at each genomic segment;}
    \item{\code{seg}}{the segment number; and}
    \item{\code{parent.index}}{the index of the clone from which this clone is descended (equals 0 if
    the clone is an original tumor clone).}
  }

  Each list in the \code{seq} component contains seven columns:
  \describe{
    \item{\code{chr}}{the chromosome number;}
    \item{\code{start}}{the locus of the simulated SNVs;}
    \item{\code{seg}}{the segment on which each SNV occurs;}
    \item{\code{mut.id}}{the id unique id number for each simulated SNV;}
    \item{\code{mutated.copies}}{the number of copies of the mutated
      allele at each SNV;}
    \item{\code{alllele}}{which allele (A or B) is mutated at each SNV; and}
    \item{\code{normal.copies}}{the number of copies of the unmutated
      allele at each SNV.}
  }

  The \code{cn.data} component contains seven columns:
  \describe{
    \item{\code{chr}}{the chromosome number;}
    \item{\code{seq}}{a unique segment identifier;}
    \item{\code{LRR}}{simulated segment-wise log ratios;}
    \item{\code{BAF}}{simulated segment-wise B allele frequencies;}
    \item{\code{X} and \code{Y}}{simulated intensities for two separate
      alleles/haplotypes per segment; and}
    \item{\code{markers}}{the simulated number of SNPS per segment.}
  }

  The \code{seq.data} component contains eight columns:
  \describe{
    \item{\code{chr}}{the chromosome number;}
    \item{\code{seq}}{a unique "segment" identifier;}
    \item{\code{mut.id}}{a unique mutation identifier;}
    \item{\code{refCounts} and \code{varCounts}}{the simulated numbers
      of reference and variant counts per mutation;}
    \item{\code{VAF}}{the simulated variant allele frequency;}
    \item{\code{totalCounts}}{the simulated total number of read counts;
      and}
    \item{\code{status}}{a character (that should probably be a factor)
      indicating whether a variant should be viewed as somatic or
      germline.}
  }

}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  Mark Zucker \email{zucker.64@buckeyemail.osu.edu}
}
\examples{
# Set of 300 simulated 'tumors' generated by artificially mixing and
# altering real data; 60 samples with one #clone, 60 with 2 clones,
# ..., 60 with 5 clones.
data("hapmapSegments", package = "CloneData")
data("snpPositions", package = "CloneData")
IDset <- c('NA07019', 'NA12234', 'NA12249', 'NA12753', 'NA12761',
           'NA18545', 'NA18975', 'NA18999', 'NA18517')
# Generating the data set:
\dontrun{
generateMixtures(dataPath = 'mixdat', mixPath = 'mixsim',
                 nPerK = rep(60,5),  segmentedData = hapmapSegments,
                 ID_pool = IDset, pos = snpPositions)
}
}
\keyword{array}
