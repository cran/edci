\name{oregMclust}
\alias{oregMclust}
\alias{regparm}
\alias{plot.oregMclust}
\alias{print.oregMclust}
\title{Orthogonal Regression Clustering}
\description{
  Computation of center points for regression data by means of orthogonal
  regression. A cluster method based on redescending M-estimators is
  used.
}

\usage{
  oregMclust(datax, datay, bw, method = "const",
    xrange = range(datax), yrange = range(datay),
    prec = 4, na = 1, sa = NULL, nl = 10, nc = NULL,
    brmaxit = 1000)

  regparm(reg)

  \method{plot}{oregMclust}(x, datax, datay, prec = 3, rcol = "black",
  rlty = 1, rlwd = 3, ...)

  \method{print}{oregMclust}(x, ...)
}

\arguments{
  \item{datax, datay}{numerical vectors of coordinates of the
    observations. Alternatively, a matrix with two or three columns can
    be given. Then, the first two columns are interpreted as coordinates of
    the observations and, if available, the third is passed to parameter
    \code{sa}.}
  \item{bw}{positive number. Bandwidth for the cluster method.}
  \item{method}{optional string. Method of choosing starting values for
    maximization. Possible values are:
    \itemize{
      \item "const": a constant number of angles for every observation
      is used. By default, one horizontal line through any observation
      is used as starting value. If a value for parameter \code{na} is passed,
      \code{na} lines through any observation are used. Alternatively, with the
      parameter \code{sa} a proper starting angle for every observation can
      be specified. In this case, \code{na} is ignored. The length of
      \code{sa} must be the number of observations.
      \item "all": every line through any two observations is used.
      \item "prob": Clusters are searched iteratively with randomly
      chosen starting values until either no new clusters are found
      (default), or until \code{nc} clusters are found. The precision of
      distinguishing the clusters can be tuned with the parameter
      \code{prec}. In each iteration, \code{nl} times a line through two
      randomly chosen observations is used as starting value.
    }
  }
  \item{xrange, yrange}{optional numerical intervals describing the
    domains of the observations. This is only used for normalization of the
    data. Note that both intervals should have approximately the same
    length or should be transformed otherwise. This is not done
    automatically, since this transformation affects the choice of the
    bandwidth.}
  \item{prec}{optional positive integer. Tuning parameter for
    distinguishing different clusters, which is passed to
    \code{\link{deldupMclust}}.}
  \item{na}{optional positive integer. Number of angles per observation
    used as starting values for \code{method = "const"} (default).} 
  \item{sa}{optional numerical vector. Angles (within \code{[0,2pi)}) used as
    starting values for \code{method = "const"} (default).}
  \item{nl}{optional positive integer. Number of starting lines in each
    iteration for \code{method = "prob"}.}
  \item{nc}{optional positive integer. Number of clusters to search if method
    \code{"const"} is chosen. Note that if \code{nc} is too
    large, i.e., \code{nc} clusters cannot be found, the function does not
    terminate. Attention! Using Windows, it is impossible to interrupt the routine
    manually in this case!}
  \item{brmaxit}{optional positive integer. Since the maximization could
    be very slow in some cases depending on the starting value, the
    maximization is stopped after \code{brmaxit} iterations.}

  \item{reg, x}{object returned from \code{oregMclust}.}

  \item{rcol, rlty, rlwd}{optional graphic parameters used for plotting
    regression lines.}
  \item{...}{additional parameters passed to \code{plot}.}
}
\details{
  \code{oregMclust} implements a cluster method based on redescending
  M-estimators for the case of orthogonal regression. This method is
  introduced by Mueller and Garlipp in 2003 (see references).

  \code{regparm} transforms the columns "alpha" and "beta" to
  "intersept" and "slope".

  See also \code{\link{bestMclust}}, \code{\link{projMclust}}, and
  \code{\link{envMclust}} for choosing the 'best' clusters out of all
  found clusters.
}

\value{
  A numerical matrix containing one row for every found
  regression center line. The columns "alpha" and "beta" are their
  parameters in the representation (cos(alpha), sin(alpha)) * (x,y)' = beta, where alpha is within \code{[0,2pi)}. For the alternative
  representation y = mx + b, the return value can be passed to \code{regparm}.

  The columns "value" and "count" give the value of the objective
  function and the number how often they are found.
}
\references{
  Mueller, C. H., & Garlipp, T. (2005). Simple consistent cluster methods based on redescending M-estimators with an application to edge identification in images. Journal of Multivariate Analysis, 92(2), 359--385.
}

\author{Tim Garlipp, \email{TimGarlipp@gmx.de}}

\seealso{\code{\link{bestMclust}}, \code{\link{projMclust}},
  \code{\link{envMclust}}, \code{\link{deldupMclust}}}

\examples{
  x = c(rnorm(100, 0, 3), rnorm(100, 5, 3))
  y = c(-2 * x[1:100] - 5, 0.5 * x[101:200] + 30)/2
  x = x + rnorm(200, 0, 0.5)
  y = y + rnorm(200, 0, 0.5)

  reg = oregMclust(x, y, 1, method = "prob")
  reg = projMclust(reg, x, y)
  reg
  plot(bestMclust(reg, 2, crit = "proj"), x, y)
}

\keyword{ cluster }
\keyword{ regression }
\keyword{ robust }
