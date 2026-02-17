\name{phf}
\alias{phf}
\title{Half-sib Family Phasing}
\description{
Phases a half-sib family using the block structure and an imputed sire haplotype matrix.
}
\usage{
phf(GenotypeMatrix, blockMatrix, sirePhasedMatrix)
}
\arguments{
  \item{GenotypeMatrix}{
    \code{matrix}. Half-sib genotypes (one half-sib per row; SNPs ordered by mapping position in the columns).
    Data should be numeric: \code{0, 1, 2} for \code{AA, AB, BB}. Use \code{9} for missing data.
  }
  \item{blockMatrix}{
    \code{matrix}. Blocking structure (output of \code{\link{bmh}}).
  }
  \item{sirePhasedMatrix}{
    \code{matrix}. Imputed sire haplotypes (output of \code{\link{ssp}}).
  }
}
\value{
Returns a matrix containing the phased parental haplotypes of the half-sibs (\strong{two rows per individual}).
Alleles are coded as \code{0} (A), \code{1} (B), and \code{9} (missing/unphased).
}
\note{
The genotype matrix must contain individuals from one half-sib family and one ordered chromosome.
This function is used by \code{\link{aio}} for complete phasing of a half-sib group.
}
\seealso{
\code{\link{aio}}
}
\examples{
genotype <- matrix(c(
  2,1,0,
  2,0,0,
  0,0,2
), byrow = TRUE, ncol = 3)

block <- bmh(genotype)
phf(genotype, block, ssp(block, genotype))
}
\keyword{snp}
\keyword{phase}
\keyword{haplotype}
\keyword{inference}
\keyword{reconstruction}
