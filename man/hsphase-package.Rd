\name{hsphase-package}
\alias{hsphase-package}
\alias{hsphase}
\docType{package}
\title{
Phasing, Pedigree Reconstruction, Sire Imputation and Recombination Events for Half-sib Families
}
\description{
Identification of recombination events, haplotype reconstruction and sire imputation using half-sib family SNP data.
}
\details{
\tabular{ll}{
Package: \tab hsphase\cr
Type: \tab Package\cr
Version: \tab 3.0.0\cr
Date: \tab 2026-02-15\cr
License: \tab GPL-3\cr
}
\strong{Main functions:}\cr
\code{\link{bmh}}: Block partitioning\cr
\code{\link{ssp}}: Sire inference\cr
\code{\link{aio}}: Phasing\cr
\code{\link{imageplot}}: Image plot of the block structure\cr
\code{\link{rpoh}}: Reconstruct pedigree based on opposing homozygotes\cr

\strong{Auxiliary functions:}\cr
\code{\link{hss}}: Half-sib family splitter\cr
\code{\link{cs}}: Chromosome splitter\cr
\code{\link{para}}: Parallel data analysis\cr
}
\author{
Mohammad H. Ferdosi \email{mferdosi@une.udu.au}, Cedric Gondro \email{cgondro2@une.edu.au}\cr
Maintainer: Mohammad H. Ferdosi \email{mhferdosi@yahoo.com}
}
\references{
Ferdosi, M. H., Kinghorn, B. P., van der Werf, J. H., & Gondro, C. (2013).
Effect of genotype and pedigree error on detection of recombination events, sire imputation and haplotype inference using the hsphase algorithm.
In \emph{Proc. Assoc. Advmt. Anim. Breed. Genet} (Vol. 20, pp. 546--549). AAABG; Napier, New Zealand.
\cr\cr
Ferdosi, M. H., Kinghorn, B. P., van der Werf, J. H. J., & Gondro, C. (2014).
Detection of recombination events, haplotype reconstruction and imputation of sires using half-sib SNP genotypes.
\emph{Genetics Selection Evolution}, 46(1), 11.
\cr\cr
Ferdosi, M. H., Kinghorn, B. P., van der Werf, J. H. J., Lee, S. H., & Gondro, C. (2014).
hsphase: an R package for pedigree reconstruction, detection of recombination events, phasing and imputation of half-sib family groups.
\emph{BMC Bioinformatics}, 15(1), 172.
\cr\cr
Ferdosi, M. H., & Boerner, V. (2014).
A fast method for evaluating opposing homozygosity in large SNP data sets.
\emph{Livestock Science}.\cr\cr
Sahoo S., Ferdosi M.H., van der Werf J.H.J., and de las Heras-Saldana S., et al. (2025) \emph{Proc. Assoc. Advmt. Anim. Breed. Genet}. 26: 323
}
\examples{
genotype <- matrix(c(
  0,0,0,0,1,2,2,2,0,0,2,0,0,0,
  2,2,2,2,1,0,0,0,2,2,2,2,2,2,
  2,2,2,2,1,2,2,2,0,0,2,2,2,2,
  2,2,2,2,0,0,0,0,2,2,2,2,2,2,
  0,0,0,0,0,2,2,2,2,2,2,0,0,0
), ncol = 14, byrow = TRUE)

ssp(bmh(genotype), genotype)
aio(genotype)
imageplot(bmh(genotype), title = "ImagePlot example")
rplot(genotype, 1:14)
}
\keyword{block}
\keyword{halfsib}
\keyword{sire}
\keyword{phasing}
