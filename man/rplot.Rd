\name{rplot}
\alias{rplot}
\title{Recombination Plot}
\description{
Creates a plot showing the sum of recombination events across a half-sib family.
}
\usage{
rplot(x, distance, start = 1, end = ncol(x), maximum = 100,
      overwrite = FALSE, method = "constant")
}
\arguments{
  \item{x}{
    \code{matrix}. Half-sib genotypes (one half-sib per row; SNPs ordered by mapping position in columns).
    Numeric coding: \code{0, 1, 2} for \code{AA, AB, BB}. Use \code{9} for missing data.
  }
  \item{distance}{
    \code{numeric} (or \code{integer}). Physical distances between markers (length must match \code{ncol(x)} or the plotted range).
  }
  \item{start}{
    \code{integer}. First marker index for the plot.
  }
  \item{end}{
    \code{integer}. Last marker index for the plot.
  }
  \item{maximum}{
    \code{integer}. Maximum number of recombinations to show (higher recombination rates will be omitted).
  }
  \item{overwrite}{
    \code{logical}. Draw over the current plot (default \code{FALSE}).
  }
  \item{method}{
    \code{character}. Method passed to \code{\link{pm}} (e.g., \code{"constant"} or \code{"relative"}).
  }
}
\examples{
genotype <- matrix(c(
  0,2,0,1,0,
  2,0,1,2,2,
  2,2,1,0,2,
  2,2,1,1,1,
  0,0,2,1,0
), ncol = 5, byrow = TRUE)

rplot(genotype, c(1,2,3,4,8))
}
\keyword{snp}
\keyword{halfsib}
\keyword{phase}
\keyword{recombination}
