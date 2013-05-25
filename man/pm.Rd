\name{pm}
\alias{pm}

\title{
Probability Matrix
}
\description{
Creates a recombination matrix based on the blocking structure.} 
\usage{
pm(blockMatrix) 
}

\arguments{
  \item{blockMatrix}{ 
\code{matrix} blocking structure (Output of \code{\link{bmh}})
}
}
\details{
This function finds the recombination between two consecutive sites, and marks the recombination site with a 1; if there are unknown sites between two blocks it will also mark these sites with a 1.}
\examples{
genotype <- matrix(c(
0,2,0,1,0,
2,0,1,2,2,
2,2,1,0,2,
2,2,1,1,1,
0,0,2,1,0), ncol=5, byrow=TRUE) # each row contains the SNP of individuals

(result <- bmh(genotype))
pm(result)
}
