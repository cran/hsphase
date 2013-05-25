\name{hsphase-package}
\alias{hsphase-package}
\alias{hsphase}
\docType{package}
\title{
Phasing, Sire Imputation and Recombination Events for Half-sib Families
}
\description{
Identification of recombination events, haplotype reconstruction and sire imputation using half-sib family SNP data.
}
\details{
\tabular{ll}{
Package: \tab hsphase\cr
Type: \tab Package\cr
Version: \tab 1.2.5\cr
Date: \tab 2013-10-30\cr
License: \tab GPL 3 \cr
}
\strong{Main Functions:} \cr
\code{\link{bmh}}: Block partitioning\cr
\code{\link{ssp}}: Sire inference\cr
\code{\link{aio}}: Phasing\cr
\code{\link{imageplot}}: Image plot of the block structure \cr
\code{\link{rpoh}}: Reconstruct pedigree based on opposing homozygote

\strong{Auxiliary Functions} \cr
\code{\link{hss}}: Half-sib family splitter\cr
\code{\link{cs}}: Chromosome splitter\cr
\code{\link{para}}: Parallel Data Analysis\cr


\emph{Note}: These functions can be used to analyse large datasets.

}
\author{
Mohammad H. Ferdosi <mferdosi@une.udu.au>, Cedric Gondro <cgondro2@une.edu.au>

Maintainer: Mohammad H. Ferdosi <mferdosi@myune.udu.au>
}
\references{
Ferdosi MH, Kinghorn BP, JHJ van der Werf, Gondro C: Effect of genotype and pedigree error on block partitioning, sire imputation and haplotype inference using the hsphase algorithm. In: AAABG; Napier, New Zealand.  2013.
\cr\cr
Ferdosi, MH, BP Kinghorn, JHJ van der Werf and C Gondro (2013). Detection of recombination events, haplotype reconstruction and imputation of sires using half-sib SNP genotypes. under review
\cr\cr
Ferdosi, MH, BP Kinghorn, JHJ van der Werf and C Gondro (2013). HSPhase: an R package for detection of recombination events, phasing and imputation of half-sib family groups. under review 
\cr\cr

}

\keyword{block,half-sib,sire inference,phase }

\examples{
genotype <- matrix(c(
  0,0,0,0,1,2,2,2,0,0,2,0,0,0,
  2,2,2,2,1,0,0,0,2,2,2,2,2,2,
  2,2,2,2,1,2,2,2,0,0,2,2,2,2,
  2,2,2,2,0,0,0,0,2,2,2,2,2,2,
  0,0,0,0,0,2,2,2,2,2,2,0,0,0), ncol=14, byrow=TRUE) 
ssp(bmh(genotype), genotype)
aio(genotype)
imageplot(bmh(genotype), title="ImagePlot example")
rplot(genotype, c(1:14))
}