\name{impute}
\alias{impute}
\title{
Impute of Low Density SNP Marker to High Density (Paternal Strand)}
\description{
Impute the paternal strand from low density to high density utilising high density sire haplotype.}
\usage{
impute(halfsib_paternal_ld, halfsib_block_ld, sire_ld, sire_hd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{halfsib_paternal_ld}{\code{Matrix}
half-sib paternal strand (output of {\code{phf}} - the matrix should have rownames which are the sample IDs and colnames which are the SNP names)
}
  \item{halfsib_block_ld}{\code{Matrix}
halfsib block structure (Output of \code{bmh})}
  \item{sire_ld}{\code{Matrix}
imputed sire (output of {\code{ssp}}  - the matrix should have rownames which are the sample IDs and colnames which are the SNP names)}
  \item{sire_hd}{\code{Matrix}
haplotype of sire (this parameter can be sequence data or any phased sire - the matrix should have rownames which are the sample IDs and colnames which are the SNP names)}
}

\value{
Return an imputed half-sib matrix. 
}

\seealso{
 \code{\link{bmh}},  \code{\link{ssp}} and  \code{\link{phf}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{impute}
\keyword{Low_Density}
\keyword{High_Density}
