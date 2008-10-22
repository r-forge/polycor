\name{polyserial}
\alias{polyserial}

\title{Polyserial Correlation}

\description{
  Computes the polyserial correlation (and its standard error)
  between a quantitative variable and an ordinal variables, based on the assumption
  that the joint distribution of the quantitative variable and a latent continuous variable underlying the ordinal 
  variable is bivariate normal. Either
  the maximum-likelihood estimator or a quicker ``two-step'' approximation is available. For the ML
  estimator the estimates of the thresholds and the covariance matrix of the estimates are also available.
}
\usage{
polyserial(x, y, ML = FALSE, control = list(), std.err = FALSE, maxcor=.9999, bins=4)
}

\arguments{
  \item{x}{a numerical variable.}
  \item{y}{an ordered categorical variable; can be numeric, a factor,
    or an ordered factor, but if a factor, its levels should be in proper order.}
  \item{ML}{if \code{TRUE}, compute the maximum-likelihood estimate; if \code{FALSE}, the default, compute a quicker
    ``two-step'' approximation.}
  \item{control}{optional arguments to be passed to the \code{optim} function.}
  \item{std.err}{if \code{TRUE}, return the estimated variance of the correlation (for the two-step estimator)
    or the estimated covariance matrix of the correlation and thresholds (for the ML estimator);
    the default is \code{FALSE}.}
  \item{maxcor}{maximum absolute correlation (to insure numerical stability).}
  \item{bins}{the number of bins into which to dissect \code{x} for a test of bivariate normality; the default is 4.}
}

\value{
  If \code{std.err} is \code{TRUE}, 
  returns an object of class \code{"polycor"} with the following components:
  \item{type}{set to \code{"polyserial"}.}
  \item{rho}{the polyserial correlation.}
  \item{cuts}{estimated thresholds for the ordinal variable (\code{y}), for the ML estimator.}
  \item{var}{the estimated variance of the correlation, or, for the ML estimator, \
    the estimated covariance matrix of the correlation and thresholds.}
  \item{n}{the number of observations on which the correlation is based.}
  \item{chisq}{chi-square test for bivariate normality.}
  \item{df}{degrees of freedom for the test of bivariate normality.}
  \item{ML}{\code{TRUE} for the ML estimate, \code{FALSE} for the two-step estimate.}
  Othewise, returns the polyserial correlation.
}

\references{
    Drasgow, F. (1986) 
    Polychoric and polyserial correlations. 
    Pp. 68--74 in S. Kotz and N. Johnson, eds.,
    \emph{The Encyclopedia of Statistics, Volume 7.} Wiley.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{hetcor}}, \code{\link{polychor}}, \code{\link{print.polycor}},
  \code{\link[stats]{optim}}}

\examples{
set.seed(12345)
data <- rmvnorm(1000, c(0, 0), matrix(c(1, .5, .5, 1), 2, 2))
x <- data[,1]
y <- data[,2]
cor(x, y)  # sample correlation
y <- cut(y, c(-Inf, -1, .5, 1.5, Inf))
polyserial(x, y)  # 2-step estimate
polyserial(x, y, ML=TRUE, std.err=TRUE) # ML estimate
}

\keyword{models}

