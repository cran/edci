\name{bestMclust}
\alias{bestMclust}
\alias{projMclust}
\alias{envMclust}
\title{Choose 'best' clusters}
\description{
  
  Chooses the 'best' regression cluster(s), if the number of true
  clusters is known.

}

\usage{

  bestMclust(clust, nc=1, crit="value")
  projMclust(clust, x, y)
  envMclust(clust, x, y, dist=0)

}

\arguments{

  \item{clust}{Cluster object returned by \code{oregMclust} or \code{circMclust}.}
  \item{nc}{Number of 'best' clusters.}
  \item{crit}{Name of the column to determine the best clusters.}
  \item{x,y}{Original observations.}
  \item{dist}{Maximal distance of observation from cluster center.}

}
\details{
  \code{oregMclust} and \code{circMclust} return a matrix containing not
  only the parameters of the found clusters but the value of the heights
  of the corresponding local maxima as well as how often each cluster is
  found. Both are reasonable criteria for choosing 'best' clusters, which
  can be done by \code{bestMclust}.
  Additional criteria could be the number of observations projected to
  each cluster or the number of observations lying in a certain
  neighbourhood of the cluster center point.

  \code{projMclust} adds a column \code{proj} to \code{clust} which
  contains the number of points belonging to each cluster in the sense
  that each observation belongs to the cluster with shortest orthogonal
  distance. If \code{clust} is comming from \code{circMclust} a second
  column \code{projrel} is added which contains this number relativ to
  the radius of the particular circle.

  \code{envMclust} adds a column \code{env} to \code{clust} which
  contains the number of observations lying in a
  \code{dist}-neighbourhood of each cluster center. If \code{clust} is
  comming from \code{circMclust} a second column \code{envrel} is added
  which contains this number relativ to the radius of the particular
  circle.
}
\value{
  Both functions return a matrix of clusters.
}
\references{
  Mueller, C.H., Garlipp, T. (2003) Simple consistent cluster methods
  based on redescending M-estimators with an application to edge
  identification in images, to appear in \emph{JMVA}.
}
\author{Tim Garlipp, \email{TimGarlipp@gmx.de}}
\keyword{ cluster }
\keyword{ regression }
\keyword{ robust }
