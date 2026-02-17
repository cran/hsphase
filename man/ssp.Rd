\name{ssp}
\alias{ssp}
\title{Sire Imputation and Phasing}
\description{
Infers (imputes) and phases the sire haplotypes based on the block structure matrix and homozygous sites of the half-sib genotype matrix.
}
\usage{
ssp(blockMatrix, genotypeMatrix)
}
\arguments{
  \item{blockMatrix}{
    \code{matrix}. Block structure (output of \code{\link{bmh}}).
  }
  \item{genotypeMatrix}{
    \code{matrix}. Half-sib genotypes (one individual per row).
    Genotypes coded as \code{0, 1, 2} for \code{AA, AB, BB}. Use \code{9} for missing data.
  }
}
\value{
A \code{matrix} with two rows (one per sire haplotype) and columns corresponding to SNPs in genotype order.
Alleles are coded as \code{0} (A) and \code{1} (B). Alleles that could not be imputed are coded as \code{9}.
}
\seealso{
\code{\link{phf}}, \code{\link{aio}}, \code{\link{imageplot}}
}
\examples{
genotype <- matrix(c(
  0,2,1,1,1,
  2,0,1,2,2,
  2,2,1,0,2,
  2,2,1,1,1,
  0,0,2,1,0
), ncol = 5, byrow = TRUE)

ssp(bmh(genotype), genotype)
}
\keyword{snp}
\keyword{halfsib}
\keyword{phase}
