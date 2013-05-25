\name{rpoh}
\alias{rpoh}
\title{
Reconstruct Pedigree Based on Matrix of Opposing Homozygotes}
\description{
Reconstructs a half-sib pedigree based on a matrix of opposing homozygotes.}
\usage{
rpoh(genotypeMatrix, oh, forwardVectorSize = 30, excludeFP = TRUE, nsap = 3,
maxRec = 15, intercept = 26.3415, coefficient = 77.3171, snpnooh, method, maxsnpnooh)
}
\arguments{
  \item{genotypeMatrix}{
\code{matrix} genotype of one chromosome (data should be numeric. Use 0, 1 and 2 for respectively AA, AB and BB. Use 9 for missing data)
}
\item{oh}{
  \code{integer} Opposing homozygotes matrix (Output of \code{\link{ohg}}) }


  \item{forwardVectorSize}{
   \code{integer} number of heterozygous sites used to validate recombination events or check for genotyping errors}


\item{excludeFP}{
  \code{logical} excludes SNPs that may cause heterozygous sites in the sire due to genotyping errors or map errors}


\item{nsap}{
  \code{integer} number of SNP per block to validate recombinations}

\item{maxRec}{
  \code{integer} maximum number of expected recombinations per individual}
  
\item{intercept}{
  \code{integer} intercept of fitted model}
  
\item{coefficient}{
  \code{integer} coefficient of fitted model}
  
  \item{snpnooh}{
  \code{Integer} number of SNPs used to create \emph{oh} matrix (this number must be divided by 1000)}
  
\item{method}{
  \code{Character} pedigree reconstruction method}
  
\item{maxsnpnooh}{
  \code{numeric} the maximum number of allowing opposinghomozygote in a half-sib family
}
}

\details{
Four methods \emph{simple}, \emph{recombinations}, \emph{calus} and \emph{manual} can be
utilised to reconstruct the pedigree. The following example shows the arguments require for each method.

pedigree <- rpoh(oh = ohg(genotype), snpNooh = 50, method = "simple")
pedigree <- rpoh(genotypeMatrix = genotype, oh = ohg(genotype), maxRec = 10 
, method = "recombinations")
pedigree <- rpoh(genotypeMatrix = genotype, oh = ohg(genotype),method = "calus")
pedigree <- rpoh(oh = ohg(genotype), maxsnpnooh = 31662, method = "manual")
}

\value{
Returns a data frame with two columns, the first column is animals' ID and the second column is sire identifiers (randomly generated).
}
\note{
Method can be \emph{recombinations} or \emph{simple}. Please refer to vignette for more information.
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