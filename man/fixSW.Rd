\name{fixSW}
\alias{fixSW}
\title{Fixing Switch Errors}
\description{
Fix switch errors in haplotypes for a half-sib family.
}
\usage{
fixSW(haplotype, ohMax = 0, windowsSize = 100, minLength = 100, cpus = 2)
}
\arguments{
  \item{haplotype}{
    \code{matrix}. Haplotypes for a half-sib family (two rows per individual).
  }
  \item{ohMax}{
    \code{integer}. Maximum tolerated opposing homozygotes when grouping each partition (increase if genotyping errors exist).
  }
  \item{windowsSize}{
    \code{integer}. Partition size (number of SNPs).
  }
  \item{minLength}{
    \code{integer}. Minimum length between switches.
  }
  \item{cpus}{
    \code{integer}. Number of CPU threads.
  }
}
\value{
A haplotype \code{matrix} with switch errors corrected.
}
\seealso{
\code{\link{groupMatSingle}} and \code{\link{addSwitch}}
}
\examples{
if (requireNamespace("RhpcBLASctl", quietly = TRUE))
    RhpcBLASctl::omp_set_num_threads(2)

haplotype <- .simulateHalfsib(7, 2500, type = "haplotype")$phased
switches <- list(500,0,0,1200,c(1000,2000),500,1200)

haplotype2 <- addSwitch(haplotype, switches, 0)

gMat <- groupMatSingle(haplotype2, 100, 2, "haplotype")
imageplot(gMat, title = "Before fixing switches")

haplotype3 <- fixSW(haplotype2, 0, 100, 100)

gMat2 <- groupMatSingle(haplotype3, 100, 2, "haplotype")
imageplot(gMat2, title = "After fixing switches")
}
\keyword{haplotype}
\keyword{switch}
\keyword{fix}
