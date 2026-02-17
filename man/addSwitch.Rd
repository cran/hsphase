\name{addSwitch}
\alias{addSwitch}
\title{Add Switches}
\description{
Add switch points to haplotypes by swapping the two haplotype rows for an individual from each switch point to the end of the chromosome.
}
\usage{
addSwitch(haplotypeMatrix, switchPoints, minLength)
}
\arguments{
  \item{haplotypeMatrix}{
    \code{matrix}. Haplotypes for a half-sib family (\strong{two rows per individual}).
  }
  \item{switchPoints}{
    \code{list} of integer/numeric vectors. Length must equal the number of individuals.
    Each element contains the switch positions (0-based/1-based depends on how they were produced; see Details).
    If there are no switches for an individual, use \code{0}.
  }
  \item{minLength}{
    \code{integer}. Minimum distance between consecutive switch points.
    \strong{Note:} in the current implementation this filter may not be enforced (depends on package version).
  }
}
\details{
\strong{Important:} Each switch point causes a swap of the two haplotype rows for that individual from the switch position to the end.

The \code{switchPoints} list must have one element per individual (not per haplotype row).
If an element is \code{0}, no switch is applied for that individual.

If you rely on \code{minLength} to ignore nearby switches, verify your installed version enforces this rule.
}
\value{
A \code{matrix} of the same dimension as \code{haplotypeMatrix} with switches applied.
}
\seealso{
\code{\link{groupMatSingle}} and \code{\link{fixSW}}
}
\examples{
haplotype <- matrix(c(0, 0, 0, 0,
                      1, 1, 1, 1,
                      0, 0, 1, 1,
                      1, 1, 0, 0,
                      1, 1, 1, 1,
                      0, 0, 0, 0), byrow = TRUE, nrow = 6)

switchPoints <- list(firstInd = c(2), secondInd = c(1, 3), lastInd = 0)
addSwitch(haplotype, switchPoints, 0)
}
\keyword{switch}
\keyword{haplotype}
