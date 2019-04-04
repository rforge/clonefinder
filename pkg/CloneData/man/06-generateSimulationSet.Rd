\name{Generating simulated tumor and data sets}
\alias{generateSimulationSet}
\title{ Generating sets of simulated tumors with SNP array and SNV data }
\description{
  Generating and saving a set of simulated tumors and data that can be used in clonal heterogeneity analysis to assess accuracy of algorithms. 
}
\usage{
generateSimulationSet(simPath, dataPath, nPerK, rounds=400, nu=0,
                      pcnv=1, norm.contam=FALSE, dataPars=NULL)
}
\arguments{
  \item{simPath}{ path to which simulated tumors will be saved.}
  \item{dataPath}{ path to which simulated SNP array and/or SNV data will be saved.}
  \item{nPerK}{a vector of integers denoting the number of tumors to generate for each possible number
  of clones, where the nth entry dictates how many n-clone tumors will be generated.}
  \item{rounds}{integer; the number of branches or total 'historical' clones generated
  in the tumor simulation.}
  \item{nu}{an integer; the average number of mutations occuring per clonal branchin event.}
  \item{pcnv}{a real number between 0 to 1; the probability of a cnv occurring at each
  clonal branching event.}
  \item{norm.contam}{a logical value; determines whether to include normal 
  contamination in simulated tumor.}
  \item{dataPars}{ a list of parameters for data generation; see details .}
}
\details{
  A set of simulation can be generated including both the simulated clonally 
  heterogeneous tumors and the data generated therefrom. The size and general 
  characteristics of the tumor set, as well as the types of data to be created 
  from it (SNP array data and/or SNV data), are determined by the input parameter
  s. The script included generates three simulated data sets, each with 300 
  simulations, one with only copy number alterations (and only SNP array data), 
  one with only single nucleotide variants (SNVs) and SNV data, and one with both.
}
\value{
  The \code{generateSimulationSet} function generates and saves two
  lists for each simulation: a simulated tumor (saved in the path
  'simpath'), with objects: \code{psi}, a vector of clonal fractions,
  and \code{clones}, which is a list of tumor clones, each of which in
  turn consists of a data frame \code{cn} and a data frame \code{seq} ;
  and a simulated data object (saved in the path 'datapath'), with
  objects: \code{cn.data} and \code{se .data}. Each component is itself
  a data frame. Note that in some cases, one of these data frames may
  have zero rows or may be returned as an \code{NA}. 

  Each list in the \code{cn} component contains seven columns:
  \describe{
    \item{\code{chr}}{the chromosome number;}
    \item{\code{start}}{the starting locus of each genomic segment;}
    \item{\code{end}}{the ending locus of each genomic segment;}
    \item{\code{A}}{the first allelic copy number at each genomic segment;}
    \item{\code{B}}{the second allelic copy number at each genomic segment;} 
    \item{\code{seg}}{the segment number; and}
    \item{\code{parent.index}}{the index of the clone from which this
      clone is descended (equals 0 if the clone is an original tumor
      clone).}
  }

  Each list in the \code{seq} component contains seven columns:
  \describe{
    \item{\code{chr}}{the chromosome number;}
    \item{\code{start}}{the locus of the simulated SNVs;}
    \item{\code{seg}}{the segment on which each SNV occurs;}
    \item{\code{mut.id}}{the id unique id number for each simulated SNV;}
    \item{\code{mutated.copies}}{the number of copies of the mutated allele at each SNV;}
    \item{\code{alllele}}{which allele (A or B) is mutated at each SNV; and}
    \item{\code{normal.copies}}{the number of copies of the unmutated allele at each SNV.}
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
\references{
  Zucker MR, Abruzzo LV, Herling CD, Barron LL, Keating MJ, Abrams ZB,
  Heerema N, Coombes KR. Inferring Clonal Heterogeneity in Cancer using
  SNP Arrays and Whole Genome Sequencing. Bioinformatics. To
  appear. doi: 10.1093/bioinformatics/btz057.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  Mark Zucker \email{zucker.64@buckeyemail.osu.edu}
}
\examples{
# Simulation set with just CNVs, 300 simulations in total, 60 with 1
#clone, 60 with 2 clones... 60 with 5 clones.
\dontrun{
generateSimulationSet(simPath = 'sims-cnv', dataPath = 'data-cnv',
    nPerK = rep(60,5), rounds = 400, nu = 0, pcnv = 1, norm.contam = FALSE)
}
}
\keyword{array}
