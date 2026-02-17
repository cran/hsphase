\name{imageplot}
\alias{imageplot}

\title{Image Plot of Blocking Structure}
\description{
Create an image plot of the blocking structure.
}

\usage{
imageplot(x, title = c(), rv = FALSE, ...)
}

\arguments{
  \item{x}{
    \code{matrix}. Blocking structure (output of \code{\link{bmh}} or \code{\link{hbp}}).
  }
  \item{title}{
    \code{character} (or \code{NULL}). Title of the image plot.
  }
  \item{rv}{
    \code{logical}. If \code{TRUE}, reverse the colour scheme.
  }
  \item{\dots}{
    Optional graphical parameters. 
  }
}

\details{
White indicates regions of unknown origin; red and blue correspond to the two sire strands.
}

\author{
This is a modified version of a function written by Chris Seidel, available at
\url{http://www.phaget4.org/R/image_matrix.html}.
}

\seealso{
\code{\link{bmh}}, \code{\link{aio}}
}

\examples{
genotype <- matrix(c(
  0,2,1,1,1,
  2,0,1,2,2,
  2,2,1,0,2,
  2,2,1,1,1,
  0,0,2,1,0
), ncol = 5, byrow = TRUE)

imageplot(bmh(genotype))
}
\keyword{block}
\keyword{image}
