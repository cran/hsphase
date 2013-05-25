\name{ohg}
\alias{ohg}
\title{
Matrix of Opposing Homozygotes 
}
\description{
Creates a matrix of opposing homozygotes from the genotype matrix.
}
\usage{
ohg(genotypeMatrix, cpus = 2)
}
\arguments{
  \item{genotypeMatrix}{
{\code{matrix} genotype (Data should be numeric. Use 0, 1 and 2 respectively for AA, AB and BB. Use 9 for missing data)
}
}
 \item{cpus}{
{\code{integer} set the number of CPUs
}
}
}


\value{
Returns a square matrix (sample X sample) with the pairwise counts of opposing homozygotes.
}



\seealso{
 \code{\link{rpoh}}
}
\examples{
genotype=matrix(c(        
  2,1,0,                  
  2,0,0,                  
  0,0,2                   
  ), byrow=TRUE, ncol=3)    
  
ohg(genotype)             

}

\keyword{opposing-homozygote}
\keyword{snp}
\keyword{genoytpe}
