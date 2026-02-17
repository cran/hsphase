\name{ohd}
\alias{ohd}
\title{Opposing Homozygote Detection}
\description{
Counts, for each animal, the number of loci where it contributes opposing homozygotes in sites that imply heterozygosity in the sire.
}
\usage{
ohd(genotypeMatrix, unique_check = FALSE, SNPs = 6000)
}
\arguments{
  \item{genotypeMatrix}{
    \code{matrix}. Half-sib genotypes (one half-sib per row, SNPs ordered by mapping position in the columns).
    Data should be numeric: \code{0, 1, 2} for \code{AA, AB, BB} and \code{9} for missing.
  }
  \item{unique_check}{
    \code{logical}. If \code{TRUE}, counts opposing homozygotes using a uniqueness rule (see Details).
  }
  \item{SNPs}{
    \code{integer}. Number of SNPs to use. (Only applicable if supported by the installed version.)
  }
}
\value{
A numeric vector with the number of heterozygous sites that each sample caused.
}
\note{
This function can be used to identify pedigree errors; i.e., outliers with unusually high values.
}
\author{
This method was suggested by Bruce Tier <btier@une.edu.au> to identify pedigree errors.
}
\examples{
genotype <- matrix(c(
  2,1,0,
  2,0,0,
  0,0,2
), byrow = TRUE, ncol = 3)

ohd(genotype)
}
\keyword{opposing homozygote}
