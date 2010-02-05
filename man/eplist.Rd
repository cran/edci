\name{eplist}
\alias{eplist}
\title{Conversion of Matrices returned by edgepoints}
\description{
  The matrices returned by \code{\link{edgepoints}} are converted into a
  list of edgepoints and a list of corresponding angles.
  This is useful for processing the results of \code{\link{edgepoints}}
  by \code{\link{oregMclust}}.
}

\usage{
  eplist(data, maxval, test=FALSE, xc=NULL, yc=NULL)
}

\arguments{
  \item{data}{list object returned from \code{\link{edgepoints}}.}
  \item{maxval}{positive numbers. Critical value for deciding, if an
    pixel belongs to an edge or not.}
  \item{test}{optional boolean. Must be set to \code{TRUE} if
    \code{edgepoints} was used with \code{estimator="test_mean"} or
    \code{estimator="test_median"}. Then \code{maxval} is the level of
    the test.}
  \item{xc,yc}{optional numerical vectors defining the coordinates of
    the edgepoints. A pixel with jump height data[[1]][i,j] gets the
    coordinates (xc[i],yc[j]). As default, the coordinates are assumed as
    equidistant within [0,1], e.g. if the matrix of jump heights is a
    n-m-matrix, the pixel at position i,j gets the coordinates
    (i/n,j/m).}
}

\value{
  A numerical matrix is returned. The first two columns contain the
  coordinates of the pixels, for which the detected jump height is
  larger than \code{maxval} (or smaller then \code{maxval}, if
  \code{test=TRUE}). The third column contains the corresponding angles.
}

\author{Tim Garlipp, \email{TimGarlipp@gmx.de}}

\seealso{\code{\link{edgepoints}}}

\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
