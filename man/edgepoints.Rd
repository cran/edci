\name{edgepoints}
\alias{edgepoints}
\title{Edge detection in noisy images}
\description{
  Detection of edgepoints by the difference of two rotated and
  asymmetric Kernel- or M-Kernel-Estimators.
}

\usage{
  edgepoints(data, h1n, h2n,
             asteps     = 4,
             estimator  = "kernel",
             kernel     = "mean",
             score      = "gauss",
             sigma      = 1,
             kernelfunc = NULL,
             margin     = FALSE)
}

\arguments{
  \item{data}{numerical matrix describing the noisy image.}
  \item{h1n, h2n}{positive numbers. Bandwidth for the kernels.}
  \item{asteps}{optional positive integer. Number of different angles used.}
  \item{estimator}{optional string. Estimator used within the
    windows. Possible values are:
    \itemize{
      \item "kernel": Kernel-Estimators. The used kernel
      function can be choosen by \code{kernel}.
      \item "M\_mean": M-Kernel-Estimators with mean as starting value.
      The used kernel function can be choosen by \code{kernel}, the
      score function can be choosen by \code{score}.
      \item "M\_median": M-Kernel-Estimators with median as starting value.
      The used kernel function can be choosen by \code{kernel}, the
      score function can be choosen by \code{score}.
      \item "median": Median, what is a special M-Kernel-Estimator.
      \item "test\_mean": Multiple Test for equal means in both windows
      for every angle.
      \item "test\_median":  Multiple Test for equal means in both windows
      for every angle.
    }
  }
  \item{kernel}{optional string. Kernel function if
    \code{estimator="kernel"}, \code{estimator="M_mean"}, or
    \code{estimator="M_median"}. Possible values are:
    \itemize{
      \item "mean": Rectangular kernel. With \code{estimator="kernel"}
      this gives an ordinary mean estimator. With
      \code{estimator="M_mean"} or \code{estimator="M_median"} this
      gives an M-Estimator.
      \item "linear": Linear kernel function. The distance of the
      observations to the common midpoint of both windows is linearly
      measured.
      \item "linear2": Linear kernel function. The distance of the
      observations to the midpoint of the window they belong to is
      linearly measured.
      \item "gauss": Density of the normal distribution with
      sd=0.5 and zero outside [-1,1]x[-1,1].
      \item "func": Arbitrary kernel function given by \code{kernelfunc}.
    }
  }
  \item{score}{optional string. Score function for M-Kernel-Estimators
    if \code{estimator="M_mean"} or
    \code{estimator="M_median"}. Possible values are:
    \itemize{
      \item "gauss": negativ density of the normal distribution. The
      deviation can be given by parameter \code{sigma}.
      \item "huber": The Huber score function is the absolute value
      (median) within an interval [-c,c] and the square function (mean)
      outside this interval. The value of c can be given by the
      parameter \code{sigma}.
    }
  }
  \item{sigma}{optional positiv number. Parameter for the score function
    \code{"gauss"} or \code{"huber"}.}
  \item{kernelfunc}{optional function taking two numbers as arguments
    and returning a positive number. Used as kernelfunction if
    \code{kernel="func"} is given. Note that the function should be zero
    outside [-1,1]x[-1,1] and that only one function is given for both
    windows. The 'lower' part of the domain, e.g. [-1,1]x[-1,0], is used
    within one window while the 'upper' part is used within the other.}
  \item{margin}{Optional value. Results near the margin are in
    general not very reasonable. Setting \code{margin=TRUE} they are
    calculated nevertheless. With \code{margin=FALSE} the returned
    matrices have the same dimension as \code{data} but the jump heights at
    the margin are set to zero. Setting \code{margin="cut"} the retuned
    matrices are cutted down by the margins. The default is
    \code{margin=FALSE}.}
}

\details{
  \code{edgepoints} implements several versions of the RDKE method,
  introduced by Qiu in 1997.
  The original method, which uses kernel estimates, is a generalized
  version, which uses M-Kernel-Estimators, and two test procedures. The
  test procedures are multiple tests for different angles for the
  hypothesis, that there are equal means (or medians) in both windows.
  All methods apply rotating and scaling in the correct order (see
  Garlipp, 2004).
}

\value{
  Return value is a list of two numerical matrices. The first matrix
  contains the maximal jump height for every pixel, if the choosen
  estimator is not a test procedure, and p-values otherwise. The second
  matrix contains the angle, which leads to the maximal jump heigth or
  minimal p-value.
}

\references{
  Garlipp, T. (2004), On Robust Jump Detection in Regression Surface
  with Applications to Image Analysis,
  \emph{Carl-von-Ossietzky-Universitaet Oldenburg, Dissertation}
  
  Qiu, P. (1997), Nonparametric Estimation of Jump Surface, \emph{The
  Indian Journal of Statistics}, 59A, No.2, 268-294.
}

\author{Tim Garlipp, \email{TimGarlipp@gmx.de}}

\seealso{\code{\link{eplist}}}

\examples{
  y <- matrix(rep(0,60*60),nrow=60)
  y[21:40,21:40]<-1
  y <- y + matrix(rnorm(60*60,0,0.2),nrow=60)
  image(y,col=gray(seq(0,1,1/255)))

  ye <- edgepoints(y,0.05,0.05,estimator="M_median",kernel="gauss")
  image(ye[[1]]>0.7,col=gray(c(1,0)))
}

\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
