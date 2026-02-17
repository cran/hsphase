\name{groupMatSingle}
\alias{groupMatSingle}
\title{Grouping a Half-sib Family}
\description{
Group the genotype or haplotype of a half-sib family into partitions using opposing homozygotes.
}
\usage{
groupMatSingle(haplotype, windowsSize, cpus = 2, input = "haplotype", oh = 0)
}
\arguments{
  \item{haplotype}{
    \code{matrix}. Haplotypes (two rows per individual) or genotypes (one row per individual) depending on \code{input}.
  }
  \item{windowsSize}{
    \code{integer}. Partition size.
  }
  \item{cpus}{
    \code{integer}. Number of CPU threads.
  }
  \item{input}{
    \code{character}. Either \code{"haplotype"} or \code{"genotype"}.
  }
  \item{oh}{
    \code{integer}. Threshold for opposing homozygotes used for grouping (increase if genotyping errors exist).
  }
}
\value{
A grouping matrix.
}
\seealso{
\code{\link{addSwitch}} and \code{\link{fixSW}}
}
\examples{
haplotype <- .simulateHalfsib(10, 5000, type = "haplotype")$phased
gMat <- groupMatSingle(haplotype, 100, 2, "haplotype")
imageplot(gMat)
}
\keyword{group}
\keyword{matrix}
