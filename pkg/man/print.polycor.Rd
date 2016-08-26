\name{print.polycor}
\alias{print.polycor}

\title{Print Method for polycor Objects}

\description{
\code{print} method for objects of class \code{polycor}, produced by
\code{polychor} and \code{polyserial}.
}

\usage{
\method{print}{polycor}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{an object of class \code{polycor}, as returned by \code{polychor} 
    or \code{polyserial}.}
  \item{digits}{number of significant digits to be printed.}
  \item{\dots}{not used.}
}

\value{
  Invisibly returns \code{x}; used for its side effect --- i.e., printing.
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
    polychor(x, y, ML=TRUE, std.err=TRUE)  # polychoric correlation, ML estimate
    }
}

\keyword{print}
\keyword{methods}
