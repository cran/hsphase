\name{readGenotype}
\alias{readGenotype}
\title{Read and Check the Genotype File}
\description{
Reads a genotype file and optionally checks it for common formatting/data issues.
}
\usage{
readGenotype(genotypePath, separatorGenotype = " ", check = TRUE)
}
\arguments{
  \item{genotypePath}{
    \code{character}. Path to the genotype file (animals in rows and SNPs in columns).
    SNPs should be coded as \code{0, 1, 2} for \code{AA, AB, BB}. Use \code{9} for missing data.
    Please refer to the vignette for more information.
  }
  \item{separatorGenotype}{
    \code{character}. Field separator used in the genotype file.
  }
  \item{check}{
    \code{logical}. If \code{TRUE}, check the genotype file for possible errors.
  }
}
\note{
Please refer to the vignette for more information.
}
\value{
A genotype \code{matrix}.
}
\keyword{io}
\keyword{genotype}
