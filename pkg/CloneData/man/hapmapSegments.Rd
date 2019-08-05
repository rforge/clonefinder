\name{hapmapSegments}
\alias{hapmapSegments-data}
\alias{hapmapSegments}
\docType{data}
\title{Segmented HapMap Copy Number Data}
\description{
  Data from 225 HapMap control samples that has been analyzed with
  DNAcopy to identify segments and record information about copy number
  variaiton.
}
\usage{
data(hapmapSegments)
}
\format{A data frame (\code{hapmapSegments}), with seven columns and 61,163
  rows. This object contains the results of performing a segmentation
  copy number analysis (using \code{DNAcopy}).
  The seven columns are:
  \describe{
    \item{\code{loc.start}}{The starting base position of the segment.}
    \item{\code{loc.end}}{The ending base postion of the segment.}
    \item{\code{seg.median}}{The median log R ratio across the segment.}
    \item{\code{SamID}}{The HapMap sample ID.}
    \item{\code{chrom}}{The chromosome on which the segment is located,
      stored as a single character in \{1, 2, .., 22, X, Y\}.}
    \item{\code{AvgBAF}}{The average B allele frequency across the
      segment.} 
    \item{\code{num.mark}}{The number of markers (i.e., measured SNPs)
      located in the segment.}
  }
}
\source{
  BeadChip readings derived from 225 HapMap controls assessed on
  Human610-Quadv1 BeadChips were downloaded from the Gene Expression
  Omnibus (http://www.ncbi.nlm.nih.gov/geo; accession number GSE17205,
  73 CEU samples; accession number GSE17206, 75 CH + JP; accession
  number GSE17207, 77 YRI). Raw BeadChip data from 168 patients with CLL
  and 225 HapMap controls were preprocessed to decode SNP/probe
  positions and generate genotype calls, log R ratio, and B-allele
  frequency (BAF) estimates using Illumina GenomeStudio, version 2010.2
  (Illumina Inc.). Further processing to produce the segmentation
  results is described in the paper by Schweighofer et al. [1]
}
\references{
  [1] Schweighofer CD, Coombes KR, Majewski T, Barron LL, Lerner S,
  Sargent RL, O'Brien S, Ferrajoli A, Wierda WG, Czerniak BA, Medeiros LJ,
  Keating MJ, Abruzzo LV. Genomic variation by whole-genome SNP mapping
  arrays predicts time-to-event outcome in patients with chronic
  lymphocytic leukemia: a comparison of CLL and HapMap genotypes. J Mol
  Diagn. 2013 Mar;15(2):196-209.

  [2] International HapMap Consortium. The International HapMap
  Project. Nature. 2003 Dec 18;426(6968):789-96.
}
\keyword{datasets}

