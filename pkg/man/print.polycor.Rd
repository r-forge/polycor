\name{print.polycor}
\alias{print.polycor}
\alias{summary.polycor}
\alias{coef.polycor}
\alias{vcov.polycor}

\title{Standard Methods for polycor Objects}

\description{
Some standard methods for objects of class \code{polycor}, produced by
\code{\link{polychor}} and \code{\link{polyserial}}, including \code{\link{print}}, \code{\link{summary}}, \code{\link{coef}}, and \code{\link{vcov}}. The \code{summary} method simply invokes the \code{print} method.
}

\usage{
\method{print}{polycor}(x, digits = max(3, getOption("digits") - 3), ...)
\method{summary}{polycor}(object, ...)
\method{coef}{polycor}(object, correlation=TRUE, thresholds=TRUE, ...)
\method{vcov}{polycor}(object, correlation=TRUE, thresholds=TRUE, ...)
}

\arguments{
  \item{x, object}{an object of class \code{polycor}, as returned by \code{polychor} 
    or \code{polyserial}.}
  \item{digits}{number of significant digits to be printed.}
  \item{correlation}{return the estimated correlation or sampling variance of the correlation.}
  \item{thresholds}{return the estimated thresholds or sampling variances/covariances of the thresholds.}
  \item{\dots}{pass arguments from \code{summary} to \code{print}; otherwise not used.}
}


\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{polychor}}, \code{\link{polyserial}}}

\examples{
if(require(mvtnorm)){
    set.seed(12345)
    data <- rmvnorm(1000, c(0, 0), matrix(c(1, .5, .5, 1), 2, 2))
    x <- data[,1]
    y <- data[,2]
    cor(x, y) # sample correlation
    }
if(require(mvtnorm)){
    x <- cut(x, c(-Inf, .75, Inf))
    y <- cut(y, c(-Inf, -1, .5, 1.5, Inf))
    print(polychor(x, y, ML=TRUE, std.err=TRUE), 
      digits=3)  # polychoric correlation, ML estimate
    }
if(require(mvtnorm)){
    coef(polychor(x, y, ML=TRUE, std.err=TRUE))
    }
if(require(mvtnorm)){
    vcov(polychor(x, y, ML=TRUE, std.err=TRUE))
    }
}

\keyword{print}
\keyword{methods}
