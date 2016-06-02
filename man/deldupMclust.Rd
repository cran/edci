\name{deldupMclust}
\alias{deldupMclust}
\title{Delete duplicate found clusters}
\description{
  Delete clusters differing only by rounding errors or having
  maximization value zero.
}

\usage{
  deldupMclust(clust, prec = NULL,
    ncol = NULL, dz = TRUE)
}

\arguments{
  \item{clust}{numerical matrix whose columns contain the parameters of
    the clusters.}
  \item{prec}{optional positive integer. Number of decimal places for
    rounding.}
  \item{ncol}{number of columns describing the clusters. See details.}
  \item{dz}{optional boolean. With \code{dz = TRUE}, those clusters
    for which the objective function has value 0 are deleted.}
}

\details{
  Since clusters found by \code{oregMclust} or \code{circMclust} often
  differ only by rounding errors, the function
  \code{deldupMclust} can be used for rounding and deleting
  duplicates. If \code{clust} has a column named "count", its
  values are summed appropriately. Otherwise such a column is added.
  
  For parameter \code{clust} the object returned from \code{oregMclust}
  or \code{circMclust} can be used. Alternatively, an arbitrary matrix
  can be given, of which the first \code{ncol} columns describe the
  clusters.
  The parameter \code{prec} is the number of decimal places for
  rounding; the default is no rounding.
  With \code{ncol}, the number of columns that describe the clusters
  can be given. This is not needed, if \code{clust} is an object
  returned from \code{oregMclust} or \code{circMclust}.
}

\value{
  An object of the same type as \code{clust}.
}

\references{
  Mueller, C. H., & Garlipp, T. (2005). Simple consistent cluster methods based on redescending M-estimators with an application to edge identification in images. Journal of Multivariate Analysis, 92(2), 359--385.
}

\author{Tim Garlipp, \email{TimGarlipp@gmx.de}}

\seealso{\code{\link{oregMclust}}, \code{\link{circMclust}}}

\keyword{ cluster }
\keyword{ regression }
\keyword{ robust }
