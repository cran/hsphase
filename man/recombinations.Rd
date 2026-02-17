\name{recombinations}
\alias{recombinations}
\title{Recombination Number}
\description{
Counts the number of recombinations for each individual based on the block structure.
}
\usage{
recombinations(blockMatrix)
}
\arguments{
  \item{blockMatrix}{
    \code{matrix}. Block structure (output of \code{\link{bmh}}).
  }
}
\value{
A numeric vector of recombination counts with length equal to the number of individuals (rows) in \code{blockMatrix}.
}
\seealso{
\code{\link{bmh}}
}
\examples{
genotype <- matrix(c(
  2,1,0,0,
  2,0,2,2,
  0,0,2,2,
  0,2,0,0
), byrow = TRUE, ncol = 4)

recombinations(bmh(genotype))
}
\keyword{recombination}
\keyword{block}
