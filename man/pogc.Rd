\name{pogc}
\alias{pogc}
\title{Parent-Offspring Group Constructor}
\description{
Assign offspring to parents based on an opposing-homozygotes (OH) matrix.
}
\usage{
pogc(oh, genotypeError)
}
\arguments{
  \item{oh}{
    \code{matrix}. Opposing homozygotes matrix (output of \code{\link{ohg}}).
  }
  \item{genotypeError}{
    \code{integer}. Number of genotyping errors allowed when interpreting the \code{oh} matrix.
  }
}
\value{
A \code{data.frame} with two columns:
\itemize{
  \item animal ID
  \item assigned parent ID
}
}
\seealso{
\code{\link{ohg}}, \code{\link{hss}}, \code{\link{rpoh}}
}
\examples{
set.seed(1)
chr <- list()
sire <- list()

for(i in 1:5)
{
  chr[[i]] <- .simulateHalfsib(numInd = 20, numSNP = 5000, recbound = 1:10)
  sire[[i]] <- ssp(bmh(chr[[i]]), chr[[i]])
  sire[[i]] <- sire[[i]][1,] + sire[[i]][2,]
  sire[[i]][sire[[i]] == 18] <- 9
}

Genotype <- do.call(rbind, chr)
rownames(Genotype) <- 6:(nrow(Genotype) + 5)
sire <- do.call(rbind, sire)
rownames(sire) <- 1:5
Genotype <- rbind(sire, Genotype)

oh <- ohg(Genotype)
pogc(oh, 5)
}
\keyword{pedigree}
\keyword{parentage}
