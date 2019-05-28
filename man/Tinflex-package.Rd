\encoding{UTF-8}

\name{Tinflex-package}
\docType{package}

\alias{Tinflex-package}
\alias{Tinflex}

\title{Tinflex -- Universal non-uniform random number generator}

\description{
  \code{Tinflex} is a universal non-uniform random number generator
  based on the acceptence-rejection method for all distributions that
  have a piecewise twice differentiable density function.
  Required input includes the log-density function of
  the target distribution and its first and second derivatives.
}

\details{
  \tabular{ll}{
    Package: \tab Tinflex       \cr
    Type:    \tab Package       \cr
    Version: \tab 1.5           \cr
    Date:    \tab 2019-05-28    \cr
    License: \tab GPL 2 or later\cr
  }

  Package \pkg{Tinflex} serves two purposes:
  \enumerate{
    \item The installed package provides a fast routine for
    sampling from any distribution that has a piecewise twice
    differentiable density function.

    \item The \R source (including comments) presents all details of the
    general sampling method which are not entirely worked out in our
    paper cited in the see references below.
  }

  Algorithm \code{Tinflex} is a universal random variate generator based
  on transformed density rejection which is a variant of the
  acceptance-rejection method. The generator first computes
  and stores hat and squeeze functions and then uses these functions
  to generate variates from the distribution of interest. Since the
  setup procedure is separated from the generation procedure, many
  samples can be drawn from the same distribution without rerunning the
  (expensive) setup.

  The algorithm requires the following data about the distribution
  (for further details see \code{\link{Tinflex.setup}}):
  \itemize{
    \item the log-density of the targent distribution;
    \item its first and second derivative;
    \item a starting partition of its domain such that each subinterval
    contains at most one inflection point of the transformed density;
    \item a transformation for the density
    (default is the logarithm transformation).
  }    
  
  The following routines are provided.
  \describe{
    \item{\code{\link{Tinflex.setup}}}{computes hat and squeeze. The
      table is then stored in a generator object of class
      \code{"Tinflex"}.}
    \item{\code{\link{Tinflex.sample}}}{draws a random sample from a
      particular generator object.}
    \item{\code{\link{print.Tinflex}}}{prints the properties a generator
      object of class \code{"Tinflex"}.}
    \item{\code{\link{plot.Tinflex}}}{plots density, hat and squeeze
      functions for a given generator object of class \code{"Tinflex"}.}
  }

  For further details see \code{\link{Tinflex.setup}}.
}

\section{Warning }{
  It is very important to note that the user is responsible for the
  correctness of the supplied arguments. Since the algorithm works (in theory)
  for all distributions with piecewise twice differentiable density
  functions, it is not possible to detect improper arguments. It is thus
  recommended that the user inspect the generator object visually by
  means of the \code{plot} method (see \code{\link{plot.Tinflex}} for
  details).

  Package \pkg{\link[rvgtest]{rvgtest}} provides a test suite for
  non-uniform random number generators.
  (Approximate distribution functions are available through
  method \code{\link[Runuran]{pinv.new}} in package
  \pkg{\link[Runuran]{Runuran}}.)
}

\note{
  Routine \code{\link{Tinflex.sample}} is implemented both as pure \R
  code (routine \code{Tinflex.sample.R}) for documenting the algorithm
  as well as C code for fast performance.
}

\author{
  Josef Leydold \email{josef.leydold@wu.ac.at}, Carsten Botts and
  Wolfgang \enc{Hörmann}{Hoermann}.
}

\references{
  C. Botts, W. \enc{Hörmann}{Hoermann}, and J. Leydold (2013),
  Transformed Density Rejection with Inflection Points,
  Statistics and Computing 23(2), 251--260,
  DOI: 10.1007/s11222-011-9306-4
  (see also 
  Research Report Series / Department of Statistics and Mathematics
  Nr. 110, Department of Statistics and Mathematics, WU Wien,
  \url{http://epub.wu.ac.at/}).
}

\seealso{
  See \code{\link{Tinflex.setup}} for further details.
  
  Package \pkg{\link[Runuran]{Runuran}} provides a set of many other
  automatic non-uniform sampling algorithms.
}

\examples{
## Bimodal density
##   f(x) = exp( -|x|^alpha + s*|x|^beta + eps*|x|^2 )
##   with alpha > beta >= 2 and s, eps > 0

alpha <- 4.2
beta <- 2.1
s <- 1
eps <- 0.1

## Log-density and its derivatives.
lpdf   <- function(x) { -abs(x)^alpha + s*abs(x)^beta + eps*abs(x)^2 }
dlpdf  <- function(x) { (sign(x) * (-alpha*abs(x)^(alpha-1)
                           + s*beta*abs(x)^(beta-1) + 2*eps*abs(x))) }
d2lpdf <- function(x) { (-alpha*(alpha-1)*abs(x)^(alpha-2)
                          + s*beta*(beta-1)*abs(x)^(beta-2) + 2*eps) }

## Parameter cT=0 (default):
##   There are two inflection points on either side of 0.
ib <- c(-Inf, 0, Inf)

## Create generator object.
gen <- Tinflex.setup(lpdf, dlpdf, d2lpdf, ib=c(-Inf,0,Inf), rho=1.1)

## Print data about generator object.
print(gen)

## Draw a random sample
Tinflex.sample(gen, n=10)

## Inspect hat and squeeze visually in original scale
plot(gen, from=-2.5, to=2.5)
## ... and in transformed (log) scale.
plot(gen, from=-2.5, to=2.5, is.trans=TRUE)
}

\keyword{package}
\keyword{distribution}
\keyword{datagen}
