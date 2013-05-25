\name{hbp}
\alias{hbp}

\title{
Haplotype Blocks of Phased Data
}
\description{
Creates a blocking structure matrix of the half-sib family based on phased data of the sire and half-sib family.
}
\usage{
hbp(PhasedGenotypeMatrix, PhasedSireGenotype)
}

\arguments{\item{PhasedGenotypeMatrix}{\code{matrix} haplotypes for a half-sib family (two rows per individual) }
           \item{PhasedSireGenotype}{\code{matrix} haplotypes of sire}
           }
\value{
Returns a matrix where 3 or 4 stands for the SNP originating in, respectively, strands 1 and 2. 0 indicates that the source strand for the SNP is unknown.
}

\note{
The input matrices must only contain individuals from one half-sib family and one ordered chromosome.
This function is used by the \code{\link{aio}} function for complete phasing of a half-sib group.

}

\seealso{
 \code{\link{aio}}, \code{\link{ssp}}
}
\examples{
sire=matrix(c(
  0,0,0,0,0,1,                  # Haplotype one of the sire
  0,1,1,1,1,0                   # Haplotype two of the sire
  ), byrow=TRUE, ncol=6)
haplotypeHalfsib=matrix(c(
  1,0,1,1,1,1,                  # Individual one, haplotype one
  0,1,0,0,0,0,                  # Individual one, haplotype two
  0,1,1,0,1,1,                  # Individual two, haplotype one
  1,0,0,1,0,0                   # Individual two, haplotype two
  ), byrow=TRUE, ncol=6)        # 0s and 1s are alelle a and b
 hbp(haplotypeHalfsib, sire)
}
\keyword{block}
\keyword{phase}
\keyword{SNP}

