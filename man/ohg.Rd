\name{ohg}
\alias{ohg}
\title{Matrix of Opposing Homozygotes}
\description{
Creates a matrix of pairwise opposing-homozygote (OH) counts from a genotype matrix.
}
\usage{
ohg(genotypeMatrix)
}
\arguments{
  \item{genotypeMatrix}{
    \code{matrix}. Genotypes (numeric): \code{0, 1, 2} for \code{AA, AB, BB} and \code{9} for missing.
  }
}
\value{
Returns a square matrix (sample \eqn{\times} sample) of pairwise counts of opposing homozygotes.
(Some versions may return this matrix inside a named list element.)
}
\note{
This function can be slow for large datasets.
}
\author{
Ferdosi, M. H., & Boerner, V. (2014). A fast method for evaluating opposing homozygosity in large SNP data sets. Livestock Science.
}
\seealso{
\code{\link{rpoh}}
}
\examples{
genotype <- matrix(c(
  2,1,0,
  2,0,0,
  0,0,2
), byrow = TRUE, ncol = 3)

ohg(genotype)
}
\keyword{opposing-homozygote}
\keyword{snp}
\keyword{genotype}
