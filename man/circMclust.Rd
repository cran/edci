\name{circMclust}
\alias{circMclust}
\alias{plot.circMclust}
\alias{print.circMclust}
\title{Circular Clustering}
\description{
  Computation of cluster center points for circular regression data. A
  cluster method based on redescending M-estimators is used.
}

\usage{
  circMclust(datax, datay, bw,
             method="const", prec=4,
             minsx=min(datax),            maxsx=max(datax), nx=10,
             minsy=min(datay),            maxsy=max(datay), ny=10,
             minsr=0.01*max(datax,datay),
             maxsr=(max(datax,datay)-min(datax,datay)), nr=10,
             nsc=5,         nc=NULL, 
             minsd=NULL,    maxsd=NULL,
             brminx=minsx,  brmaxx=maxsx,
             brminy=minsy,  brmaxy=maxsy,
             brminr=minsr,  brmaxr=maxsr,
             brmaxit=1000)

\method{plot}{circMclust}(x, datax, datay,
                  ccol="black", clty=1, clwd=3, ...)
\method{print}{circMclust}(x, ...)
}

\arguments{
  \item{datax, datay}{numerical vectors of coordinates of the observations.}
  \item{bw}{positive number. Bandwidth for the cluster method.}
  \item{method}{optional string. Method of choosing starting values for
    maximization. Possible values are:
    \itemize{
      \item "const": a constant number of circles is used. By default,
      \code{nx}*\code{ny} equidistant midpoints within the range of the
      observations with \code{nr} different radiuses are uses as
      starting circles. The domain of the midpoints and radiuses can
      optionally be given by [\code{minsx},\code{maxsx}],
      [\code{minsy},\code{maxsy}], and [\code{minsr},\code{maxsr}].
      \item "all": every circle through any three observations is used.
      \item "prob": Clusters are searched iteratively with randomly
      chosen starting circles until either no new clusters are found
      (default), or until \code{nc} clusters are found. The precision of
      distinguishing the clusters can be tuned with the parameter
      \code{prec}. In each iteration \code{nc} times a circle through
      three randomly chosen observations is used as starting
      value. With the parameters \code{minsd} and \code{maxsd} the
      minimal and maximal distance of these observations could be
      limited.
    }
  }
  \item{nx, ny}{optional positive integer. Number of starting midpoints,
    if \code{method="const"} is chosen.}
  \item{nr}{optional positive integer. Number of starting radiuses, if
    \code{method="const"} is chosen.}
  \item{prec}{optional positive integer. Tuning parameter for
    distinguishing different clusters, which is passed to
    \code{\link{deldupMclust}}.}
  \item{minsx, maxsx, minsy, maxsy, minsr}{optional numbers
    determining the domain of starting midpoints and the range of
    radiuses, if \code{method="const"} is chosen.}
  \item{maxsr}{optional number determining the maximum radius used as
    starting value. Note, that this is valid for all \code{method}s
    while \code{minsx}, \code{maxsx}, \code{minsy}, \code{maxsy}, and
    \code{minsr} are only used by \code{method="const"}.}
  \item{nsc}{optional positive integer. Number of starting circles in each
    iteration, if \code{method="prob"} is chosen.}
  \item{nc}{optional positive integer. Number of clusters to search, if
    \code{method="prob"} is chosen. Note that, if \code{nc} is to
    large, e.g. \code{nc} clusters cannot be found, the function does not
    terminate. Attention! Using Windows, the routine cannot even be
    breaked manually in this case!}
  \item{minsd, maxsd}{optional positive numbers. Minimal and maximal
    distance of starting points which are used, if \code{method="prob"}
    is chosen.}
  \item{brminx, brmaxx, brminy, brmaxy, brminr, brmaxr}{optional
    numbers. The maximization is stopped, if the midpoint leaves the
    domain [\code{brminx},\code{brmaxx}]x[\code{brminy},\code{brmaxy}] or
    if the radius leaves [\code{brminr},\code{brmaxr}].}
  \item{brmaxit}{optional positive integer. Since the maximization could
    be very slow in some cases depending on the starting value, the
    maximization is stopped after \code{brmaxit} iterations.}
  
  \item{x}{object returned by \code{circMclust}}
  
  \item{ccol,clty,clwd}{optional graphic parameters used for plotting
    the circles.}
  \item{...}{Additional parameters passed to \code{plot}.}
}

\details{
  \code{circMclust} implements a cluster method using local
  maxima  of redescending M-estimators for the case of circular
  regression. This method is based on a method introduced by Mueller and
  Garlipp in 2003 (see references).
  
  See also \code{\link{bestMclust}}, \code{\link{projMclust}}, and
  \code{\link{envMclust}} for choosing the 'real' clusters out of the
  found.
}

\value{
  Return value is a numerical matrix containing one row for every found
  cluster circle. The columns "cx" and "cy" are their midpoints and "r"
  are the radiuses.
  
  The columns "value" and "count" give the value of the objective
  function and the number, how often they are found.
}

\references{
  Mueller, C.H., Garlipp, T. (2003) Simple consistent cluster methods
  based on redescending M-estimators with an application to edge
  identification in images, to appear in \emph{JMVA}.
}

\author{Tim Garlipp, \email{TimGarlipp@gmx.de}}

\seealso{\code{\link{bestMclust}}, \code{\link{projMclust}},
  \code{\link{envMclust}}, \code{\link{deldupMclust}}}

\examples{
  z <- (1:100 * pi)/50
  x <- c(sin(z) * 10 + 20,sin(z) * 30 + 80)+rnorm(200,0,2)
  y <- c(cos(z) * 10 + 20,cos(z) * 30 + 80)+rnorm(200,0,2)

  circ <- circMclust(x,y,5,method="prob",prec=1,nsc=20,minsd=10,maxsd=40)
  bestMclust(circ,2)
  plot(bestMclust(circ,2),x,y)
}

\keyword{ cluster }
\keyword{ regression }
\keyword{ robust }
