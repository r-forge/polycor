\name{hetcor}
\alias{hetcor}
\alias{hetcor.default}
\alias{hetcor.data.frame}
\alias{print.hetcor}
\alias{as.matrix.hetcor}
\alias{detectCores}

\title{Heterogeneous Correlation Matrix}
\description{
  \code{hetcor} computes a heterogenous correlation matrix, consisting of Pearson product-moment
  correlations between numeric variables, polyserial correlations between numeric 
  and ordinal variables, and polychoric correlations between ordinal variables.
  
  The \code{\link[parallel]{detectCores}} function is imported from the \pkg{parallel} package and re-exported.
}

\usage{
hetcor(data, ..., ML = FALSE, std.err = TRUE, 
  use=c("complete.obs", "pairwise.complete.obs"), 
  bins=4, pd=TRUE, parallel=FALSE, ncores=detectCores(logical=FALSE),
  thresholds=FALSE)
\method{hetcor}{data.frame}(data, ML = FALSE, std.err = TRUE, 
  use = c("complete.obs", "pairwise.complete.obs"), 
  bins=4, pd=TRUE, parallel=FALSE, ncores=detectCores(logical=FALSE), 
  thresholds=FALSE, ...)
\method{hetcor}{default}(data, ..., ML = FALSE, std.err = TRUE, 
  use=c("complete.obs", "pairwise.complete.obs"), 
  bins=4, pd=TRUE, parallel=FALSE, ncores=detectCores(logical=FALSE),
  thresholds=FALSE)
\method{print}{hetcor}(x, digits = max(3, getOption("digits") - 3), ...)
\method{as.matrix}{hetcor}(x, ...)
detectCores(all.tests=FALSE, logical=TRUE) 
}

\arguments{
  \item{data}{a data frame consisting of factors, ordered factors, logical variables, 
  	character variables, and/or numeric variables, or the first of several variables.}
  \item{\dots}{variables and/or arguments to be passed down.}
  \item{ML}{if \code{TRUE}, compute maximum-likelihood estimates;
    if \code{FALSE}, compute quick two-step estimates.}
  \item{std.err}{if \code{TRUE}, compute standard errors.}
  \item{bins}{number of bins to use for continuous variables in testing bivariate normality; the default is 4.}
  \item{pd}{if \code{TRUE} and if the correlation matrix is not positive-definite,
    an attempt will be made to adjust it to a
    positive-definite matrix, using the \code{\link[Matrix]{nearPD}} function in the \code{Matrix} package.
    Note that default arguments to \code{nearPD} are used (except \code{corr=TRUE}); for more control call \code{nearPD} directly.}
  \item{parallel}{if \code{TRUE} (the default is \code{FALSE}), perform parallel computations on a computer with multiple CPUs/cores.}
  \item{ncores}{the number of cores to use for parallel computations; the default is the number of physical cores detected.}
  \item{use}{if \code{"complete.obs"}, remove observations with any missing data; 
    if \code{"pairwise.complete.obs"}, compute each correlation using all observations with
    valid data for that pair of variables.}
  \item{thresholds}{if \code{TRUE} (the default is \code{FALSE}), include the estimated thresholds
  for polyserial and polychoric correlation in the returned object.}
  \item{x}{an object of class \code{"hetcor"} to be printed, or from which to extract the correlation matrix.}
  \item{digits}{number of significant digits.}
  \item{all.tests}{logical, apply all known tests; default is \code{FALSE}.}
  \item{logical}{if \code{TRUE}, detect logical CPUs/cores; if \code{FALSE}, detect physical CPUs/cores.}
}

\value{
  \code{hetcor} returns an object of class \code{"hetcor"} with the following components:
  \item{correlations}{the correlation matrix.}
  \item{type}{the type of each correlation: \code{"Pearson"}, \code{"Polychoric"}, or \code{"Polyserial"}.}
  \item{std.errors}{the standard errors of the correlations, if requested.}
  \item{n}{the number (or numbers) of observations on which the correlations are based.}
  \item{tests}{p-values for tests of bivariate normality for each pair of variables.}
  \item{NA.method}{the method by which any missing data were handled: \code{"complete.obs"} 
    or \code{"pairwise.complete.obs"}.}
  \item{ML}{\code{TRUE} for ML estimates, \code{FALSE} for two-step estimates.}
  \item{thresholds}{optionally, according to the \code{thresholds} argument, a matrix of mode list with a list of thresholds for each polyserial and polychoric correlation in the elements below the diagonal and the type of each correlation (Pearson, polyserial, or polychoric) above the diagonal.}
}

\author{John Fox \email{jfox@mcmaster.ca}}

\note{Although the function reports standard errors for product-moment correlations, transformations (the most well known
  is Fisher's \emph{z}-transformation) are available that make the approach to asymptotic normality much more rapid.}
  
\section{Warning}{
  Be careful with character variables (as opposed to factors), the values of which are
  ordered alphabetically. Thus, e.g., the values \code{"disagree"}, \code{"neutral"}, 
  \code{"agree"} are ordered \code{"agree"}, \code{"disagree"}, \code{"neutral"}.
}

\references{
    Drasgow, F. (1986) 
    Polychoric and polyserial correlations. 
    Pp. 68-74 in S. Kotz and N. Johnson, eds.,
    \emph{The Encyclopedia of Statistics, Volume 7.} Wiley.
    
    Olsson, U. (1979)
    Maximum likelihood estimation of the polychoric correlation coefficient.
    \emph{Psychometrika} \bold{44}, 443-460.
    
    Rodriguez, R.N. (1982)
    Correlation.
    Pp. 193-204 in S. Kotz and N. Johnson, eds.,
    \emph{The Encyclopedia of Statistics, Volume 2.} Wiley.
    
    Ghosh, B.K. (1966)
    Asymptotic expansion for the moments of the distribution of correlation coefficient.
    \emph{Biometrika} \bold{53}, 258-262.
    
    Olkin, I., and Pratt, J.W. (1958)
    Unbiased estimation of certain correlation coefficients.
    \emph{Annals of Mathematical Statistics} \bold{29}, 201-211.
}

\seealso{\code{\link{polychor}}, \code{\link{polyserial}}, \code{\link[Matrix]{nearPD}},
\code{\link[parallel]{detectCores}}}

\examples{
if(require(mvtnorm)){
    set.seed(12345)
    R <- matrix(0, 4, 4)
    R[upper.tri(R)] <- runif(6)
    diag(R) <- 1
    R <- cov2cor(t(R) \%*\% R)
    round(R, 4)  # population correlations
    data <- rmvnorm(1000, rep(0, 4), R)
    round(cor(data), 4)   # sample correlations
    }
if(require(mvtnorm)){
    x1 <- data[,1]
    x2 <- data[,2]
    y1 <- cut(data[,3], c(-Inf, .75, Inf))
    y2 <- cut(data[,4], c(-Inf, -1, .5, 1.5, Inf))
    data <- data.frame(x1, x2, y1, y2)
    hetcor(data)  # Pearson, polychoric, and polyserial correlations, 2-step est.
    }
if(require(mvtnorm)){
    hetcor(x1, x2, y1, y2, ML=TRUE) # Pearson, polychoric, polyserial correlations, ML est.
    }

\dontrun{
    hc <- hetcor(data, ML=TRUE)
      # parallel computation:
    hc.m <- hetcor(data, ML=TRUE, parallel=TRUE,
                   ncores=min(2, detectCores()))
    hc.m
    all.equal(hc, hc.m)
    
      # error handling:
    data$y1[data$y2 == "(0.5,1.5]"] <- NA
    hetcor(data)
    }
}

\keyword{models}

