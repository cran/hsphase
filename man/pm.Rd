\name{pm}
\alias{pm}
\title{Probability Matrix}
\description{
Creates a recombination (probability) matrix based on the blocking structure.
}
\usage{
pm(blockMatrix, method = "constant")
}
\arguments{
  \item{blockMatrix}{
    \code{matrix}. Blocking structure (output of \code{\link{bmh}}).
  }
  \item{method}{
    \code{character}. Method for creating the recombination matrix. Typically \code{"constant"} or \code{"relative"}.
  }
}
\details{
This function identifies recombination between two consecutive sites and marks recombination sites with \code{1}.
If there are unknown sites between two blocks, it marks these sites with:
\itemize{
  \item \code{1} for the \code{"constant"} method, or
  \item \eqn{1 / m} for the \code{"relative"} method, where \eqn{m} is the number of unknown sites.
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

block <- bmh(genotype)
pm(block)
}
\keyword{recombination}
\keyword{block}
