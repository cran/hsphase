\name{co}
\alias{co}
\title{Crossover Detection}
\description{
Detects possible crossover segments by comparing pairs of individuals in a half-sib family.
}
\usage{
co(genotypeMatrix)
}
\arguments{
  \item{genotypeMatrix}{
    \code{matrix}. Half-sib genotypes (one individual per row; SNPs in columns ordered by map position).
    Genotypes must be numeric: \code{0, 1, 2} for \code{AA, AB, BB} and \code{9} for missing.
  }
}
\value{
A \code{matrix} Returns a matrix with the number of crossover events for each site.
}
\examples{
genotype <- matrix(c(
  2,1,0,
  2,0,2,
  0,0,2
), byrow = TRUE, ncol = 3)

co(genotype)
}
\keyword{crossover}
