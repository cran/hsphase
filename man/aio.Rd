\name{aio}
\alias{aio}
\title{All-in-one Phasing}
\description{
Phasing of a single half-sib family group (single ordered chromosome).
}
\usage{
aio(genotypeMatrix, bmh_forwardVectorSize = 30, bmh_excludeFP = TRUE,
    bmh_nsap = 3, bmh_fillMissing = FALSE, output = "phase")
}
\arguments{
  \item{genotypeMatrix}{
    \code{matrix}. Half-sib genotypes (one half-sib per row, SNPs ordered by map position in columns).
    Data must be numeric: \code{0, 1, 2} for \code{AA, AB, BB}. Use \code{9} for missing.
  }
  \item{bmh_forwardVectorSize}{
    \code{integer}. Number of heterozygous sites used to validate recombination events or detect genotyping/map errors.
  }
  \item{bmh_excludeFP}{
    \code{logical}. Exclude SNPs that may induce false heterozygous sites in the sire due to genotyping or map errors.
  }
  \item{bmh_nsap}{
    \code{integer}. Number of SNP per block to validate recombinations (e.g. 50K \eqn{\to} 3, 700K \eqn{\to} 10).
  }
  \item{bmh_fillMissing}{
    \code{logical}. If \code{TRUE}, recombination points are placed at the mid-point of the ambiguous interval rather than marking surrounding SNPs as missing.
  }
  \item{output}{
    \code{character}. If \code{"phase"}, the function returns only the phased haplotypes matrix.
  }
}
\details{
This function calls \code{\link{bmh}}, \code{\link{ssp}}, and \code{\link{phf}}.
}
\value{
If \code{output = "phase"}, returns a haplotype matrix with \strong{two rows per individual} (first paternal, second maternal),
coded as \code{0} (allele A), \code{1} (allele B), and \code{9} (missing/unphased).

Otherwise returns a list with elements:
\itemize{
  \item \code{phasedHalfsibs}
  \item \code{sireHaplotype}
  \item \code{blockStructure}
}
}
\note{
Only this function needs to be called to phase a half-sib family.
The genotype matrix must contain individuals from a single family and a single ordered chromosome.
}
\seealso{
\code{\link{bmh}}, \code{\link{ssp}}, \code{\link{phf}}
}
\examples{
genotype <- matrix(c(       # Define a Half-sib Genotype Matrix
  2,1,0,	            # Individual 1
  2,0,0,                    # Individual 2
  0,0,2                     # Individual 3
), byrow = TRUE, ncol = 3)  # There are 3 individulas with three SNPs

aio(genotype)               # The genotypes must include only one half-sib family and one chromosome
}
\keyword{snp}
\keyword{phase}
\keyword{haplotype}
\keyword{inference}
\keyword{reconstruction}
