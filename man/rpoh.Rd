\name{rpoh}
\alias{rpoh}
\title{
Reconstruct Pedigree Based on Matrix of Opposing Homozygotes}
\description{
Reconstructs a half-sib pedigree based on a matrix of opposing homozygotes}
\usage{
rpoh(genotypeMatrix, oh, forwardVectorSize = 30,excludeFP = TRUE, nsap = 3,
maxRec = 15, intercept = 26.3415, coefficient = 77.3171, snpNooh, method)
}
\arguments{
  \item{genotypeMatrix}{
\code{Matrix} genotype of one chromosome (data should be numeric. Use 0, 1 and 2 for respectively AA, AB and BB. Use 9 for missing data).


}
\item{oh}{
  \code{Integer} Opposing homozygotes matrix (Output of \code{\link{ohg}}) }


  \item{forwardVectorSize}{
   \code{Integer} number of heterozygous sites used to validate recombination events or check for genotyping errors}


\item{excludeFP}{
  \code{logical} excludes SNPs that may cause heterozygous sites in the sire due to genotyping errors or map errors}


\item{nsap}{
  \code{Integer} number of SNP per block to validate recombinations}

\item{maxRec}{
  \code{Integer} maximum number of expected recombinations per individual}
\item{intercept}{
  \code{Integer} Intercept of fitted model}
\item{coefficient}{
  \code{Integer} coefficient of fitted model}
  \item{snpNooh}{
  \code{Integer} number of SNPs used to create \emph{oh} matrix (this number must be divided by 1000)}
\item{method}{
  \code{Character} pedigree reconstruction method}
}

\details{
This function calls the \code{\link{bmh}} and \code{\link{recombinations}} functions to count the number of recombinations in each half-sib group.   
}
\value{
Returns a data frame with two columns, the first column is animals' ID and the second column are sire identifiers (randomly generated).
}
\note{
Method can be \emph{recombinations} or \emph{simple} please refer to vignette for more information.
}



\seealso{
 \code{\link{bmh}} and \code{\link{recombinations}} 
}
\examples{
# Please run demo(hsphase)
}

\keyword{snp}
\keyword{opposing-homozygote}
\keyword{pedigree}
\keyword{reconstruction}