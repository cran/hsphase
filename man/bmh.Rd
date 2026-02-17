\name{bmh}
\alias{bmh}
\title{Block Partitioning}
\description{
Identifies the block structure (chromosome segments) in a half-sib family that each individual inherited from its sire.
}
\usage{
bmh(
  GenotypeMatrix,
  forwardVectorSize = 30,
  excludeFP = TRUE,
  nsap = 3,
  fillMissing = FALSE
)
}
\arguments{
  \item{GenotypeMatrix}{
    \code{matrix}. Half-sib genotypes (one half-sib per row; SNPs ordered by mapping position in the columns).
    Data should be numeric: \code{0, 1, 2} for \code{AA, AB, BB}. Use \code{9} for missing data.
  }
  \item{forwardVectorSize}{
    \code{integer}. Number of heterozygous sites used to validate recombination events or check for genotyping/map errors
    (50K \eqn{\to} 30, 700K \eqn{\to} 120).
  }
  \item{excludeFP}{
    \code{logical}. Exclude SNPs that may cause heterozygous sites in the sire due to genotyping errors or map errors.
  }
  \item{nsap}{
    \code{integer}. Number of SNP per block to validate recombinations (50K \eqn{\to} 3, 700K \eqn{\to} 10).
  }
  \item{fillMissing}{
    \code{logical}. Because the exact point of recombination is unknown, markers around recombination points may be set to missing.
    If \code{TRUE}, the recombination point is assumed to be in the middle of the ambiguous region, reducing missing markers.
  }
}
\value{
A matrix of block structure containing \code{1}, \code{2}, and \code{0}.
Values \code{1} and \code{2} represent the two sire strands (arbitrary labeling within a chromosome),
and \code{0} indicates unknown origin.
}
\note{
The genotype matrix must contain individuals from only one half-sib family and one ordered chromosome.
}
\seealso{
\code{\link{ssp}}, \code{\link{phf}}, \code{\link{aio}}, \code{\link{imageplot}}
}
\examples{
genotype <- matrix(c(
  0,2,1,1,1,
  2,0,1,2,2,
  2,2,1,0,2,
  2,2,1,1,1,
  0,0,2,1,0
), ncol = 5, byrow = TRUE)

bmh(genotype)
}
\keyword{snp}
\keyword{halfsib}
\keyword{phase}
