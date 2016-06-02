\name{edgecluster}
\alias{edgecluster}
\title{Edge detection in noisy images}
\description{
  \code{edgecluster} is a simple combination of \code{\link{edgepoints}}
  and \code{\link{oregMclust}}. It just passes the results of
  \code{\link{edgepoints}} to \code{\link{oregMclust}}.
}

\usage{
  edgecluster(data, h1n, h2n, maxval,
    bw = max(h1n, h2n)/qnorm(0.975),
    asteps = 4, estimator = "M_median",
    kernel = "gauss", score = "gauss",
    sigma = 1, kernelfunc = NULL)
}

\arguments{
  \item{data}{See description of \code{\link{edgepoints}}.}
  \item{h1n, h2n}{See description of \code{\link{edgepoints}}.}
  \item{asteps}{See description of \code{\link{edgepoints}}.}
  \item{estimator}{See description of \code{\link{edgepoints}}.}
  \item{kernel}{See description of \code{\link{edgepoints}}.}
  \item{score}{See description of \code{\link{edgepoints}}.}
  \item{sigma}{See description of \code{\link{edgepoints}}.}
  \item{kernelfunc}{See description of \code{\link{edgepoints}}.}
  \item{maxval}{See description of \code{\link{eplist}}.}
  \item{bw}{See description of \code{\link{oregMclust}}.}
  % The default gives the observations within the
  %'\code{edgecluster}-window' 95\% of the weight. (???)
}

\value{
  A list of two numerical matrices. The first matrix
  contains the results of \code{\link{oregMclust}}, which are the
  'edgeclusters'. The second matrix contains the result of
  \code{\link{edgepoints}}.
}

\author{Tim Garlipp, \email{TimGarlipp@gmx.de}}

\seealso{\code{\link{edgepoints}}, \code{\link{oregMclust}},
  \code{\link{eplist}}}

\examples{
  # generate a 60x60 zero matrix
  y = matrix(rep(0, 60 * 60), nrow = 60)
  # set a square-shaped set of elements to 1
  y[21:40, 21:40] = 1
  # add some noise
  y = y + matrix(rnorm(60 * 60, 0, 0.2), nrow = 60)
  # plot it
  image(y, col = gray(seq(0, 1, 1/255)))

  # find edge points of the square-shaped object 
  ec = edgecluster(y, 0.05, 0.05, 0.7,
    estimator = "M_median", kernel = "gauss")
  plot(bestMclust(ec[[1]], 4), ec[[2]], xlim = c(0, 1), ylim = c(0, 1))
}

\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
