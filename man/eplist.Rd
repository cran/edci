\name{eplist}
\alias{eplist}
\title{Conversion of matrices returned by edgepoints}
\description{
  The matrices returned by \code{\link{edgepoints}} are converted into a
  list of edge points and a list of corresponding angles.
  This is useful for processing the results of \code{\link{edgepoints}}
  by \code{\link{oregMclust}}.
}

\usage{
  eplist(data, maxval, test = FALSE, xc = NULL, yc = NULL)
}

\arguments{
  \item{data}{list object returned from \code{\link{edgepoints}}.}
  \item{maxval}{positive numbers. Critical value for deciding whether a
    pixel belongs to an edge or not.}
  \item{test}{optional boolean. Must be set to \code{TRUE} if
    \code{edgepoints} was used with \code{estimator = "test_mean"} or
    \code{estimator = "test_median"}. Then, \code{maxval} is the level of
    the test.}
  \item{xc, yc}{optional numerical vectors defining the coordinates of
    the edge points. A pixel with jump height \code{data[[1]][i,j]} gets the
    coordinates \code{(xc[i],yc[j])}. By default, the coordinates are assumed as
    equidistant within \code{[0,1]}, e.g., for an \code{(n x m)} matrix of jump heights,
    the pixel at position \code{(i,j)} gets the coordinates \code{(i/n,j/m)}.}
}

\value{
  A numerical matrix. The first two columns contain the
  coordinates of the pixels for which the detected jump height is
  larger than \code{maxval} (or smaller than \code{maxval} if
  \code{test = TRUE}). The third column contains the corresponding angles.
}

\author{Tim Garlipp, \email{TimGarlipp@gmx.de}}

\seealso{\code{\link{edgepoints}}}

\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
