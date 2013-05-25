\name{impute}
\alias{impute}
\title{
Impute of Low Density SNP Marker to High Density (Parental strand)}
\description{
Impute the parental strand from low density to high density utilising high density sire haplotype}
\usage{
impute(halfsib_parental_ld, halfsib_block_ld, sire_ld, sire_hd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{halfsib_parental_ld}{
half-sib parental strand (output of {\code{phf}} - the matrix should have rownames which are the sample IDs and colnames which are the SNP names)
}
  \item{halfsib_block_ld}{
halfsib block structure (Output of\code{bmh})}
  \item{sire_ld}{
imputed sire (output of {\code{ssp}}  - the matrix should have rownames which are the sample IDs and colnames which are the SNP names))}
  \item{sire_hd}{
haplotype of sire (this parameter can be sequence data or any phased sire  - the matrix should have rownames which are the sample IDs and colnames which are the SNP names))}
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
