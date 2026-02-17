\name{hbp}
\alias{hbp}
\title{Haplotype Blocks of Phased Data}
\description{
Creates a block-structure matrix for a half-sib family based on phased data of the sire and the half-sib family.
}
\usage{
hbp(PhasedGenotypeMatrix, PhasedSireGenotype, strand = "auto")
}
\arguments{
  \item{PhasedGenotypeMatrix}{
    \code{matrix}. Haplotypes for a half-sib family (\strong{two rows per individual}).
    Alleles should be coded as \code{0} and \code{1}; use \code{9} for missing/unphased if present.
  }
  \item{PhasedSireGenotype}{
    \code{matrix}. Haplotypes of the sire (two rows; same SNP order as \code{PhasedGenotypeMatrix}).
  }
  \item{strand}{
    \code{character}. Method for identification of paternal strand.
    Use \code{"auto"} (recommended; default) or specify \code{"1"} or \code{"2"} to force a strand definition.
  }
}
\value{
A \code{matrix} in which \code{3} or \code{4} indicates the SNP originates from, respectively, sire strand 1 or strand 2.
\code{0} indicates the origin is unknown.
}
\note{
The input matrices must contain individuals from a single half-sib family and a single ordered chromosome.
The SNP order must match between inputs.
}
\seealso{
\code{\link{aio}}, \code{\link{ssp}}
}
\examples{
sire <- matrix(c(
  0,0,0,0,0,1,                  # Haplotype one of the sire
  0,1,1,1,1,0                   # Haplotype two of the sire
  ), byrow = TRUE, ncol = 6)
  
haplotypeHalfsib <- matrix(c(
  1,0,1,1,1,1,                  # Individual one, haplotype one
  0,1,0,0,0,0,                  # Individual one, haplotype two
  0,1,1,0,1,1,                  # Individual two, haplotype one
  1,0,0,1,0,0                   # Individual two, haplotype two
  ), byrow = TRUE, ncol = 6)    # 0s and 1s are alelle a and b
  
 hbp(haplotypeHalfsib, sire)
}
\keyword{block}
\keyword{phase}
\keyword{SNP}
