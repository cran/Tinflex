\encoding{UTF-8}

\name{Tinflex.setup}
\alias{Tinflex.setup}
\alias{Tinflex.setup.C}

\title{
  Create Tinflex Generator Objects 
}

\description{
  Create a generator object of class \code{"Tinflex"} or
  \code{"TinflexC"}.
}

\usage{
Tinflex.setup(lpdf, dlpdf, d2lpdf=NULL, ib, cT=0, rho=1.1, max.intervals=1001)
Tinflex.setup.C(lpdf, dlpdf, d2lpdf=NULL, ib, cT=0, rho=1.1, max.intervals=1001)
}

\arguments{
  \item{lpdf}{log-density of targent distribution. (function)}
  \item{dlpdf}{first derivative of log-density. (function)}
  \item{d2lpdf}{second derivative of log-density. (function, optional)}
  \item{ib}{break points for partition of domain of log-density.
    (numeric vector of length greater than 1)}
  \item{cT}{parameter for transformation
    \ifelse{html}{\out{<i>T<sub>c</sub></i>}}{\eqn{T_c}}.
    (numeric vector of length 1 or of length \code{length(ib)-1})}
  \item{rho}{performance parameter:
    requested upper bound for ratio of area below hat to area
    below squeeze. (numeric)}
  \item{max.intervals}{maximal numbers of intervals. (numeric)}
}
 
\details{
  Algorithm \code{Tinflex} is a flexible algorithm that works (in
  theory) for all distributions that have a piecewise twice
  differentiable density function.
  The algorithm is based on the transformed density rejection algorithm
  which is a variant of the acceptance-rejection algorithm where
  the density of the targent distribution is transformed by means of
  some transformation
  \ifelse{html}{\out{<i>T<sub>c</sub></i>}}{\eqn{T_c}}.
  Hat and squeeze functions of the density are then constructed by means
  of tangents and secants.

  The algorithm uses family
  \ifelse{html}{\out{<i>T<sub>c</sub></i>}}{\eqn{T_c}}
  of transformations, where

  \ifelse{html}{\out{
      <table>
      <tr>
      <td><i>T<sub>0</sub> = </i>log<i>(x)</i>
      <td>for <i>c=0</i></tr>
      <tr>
      <td><i>T<sub>c</sub> = </i>sign<i>(x) x<sup>c</sup></i>
      &nbsp;&nbsp;&nbsp;&nbsp;
      <td>for <i>c!=0</i> </tr>
      </table>}}{%
    \deqn{T_c(x) = \left\{\begin{array}{l@{\quad}l}%
      \log(x) & \mbox{for $c=0\,,$}\\
      \mbox{sign}(c)\; x^c & \mbox{for $c\not=0\,.$}
      \end{array}\right.}{%
      T_0(x) = log(x)       for c=0,
      T_c(x) = sign(c) x^c  for c!=0.}
  }
  
  Parameter \eqn{c} is given by argument \code{cT}.
  
  The algorithm requires the following input from the user:
  \itemize{
    \item the log-density of the targent distribution, \code{lpdf};
    \item its first derivative \code{dlpdf};
    \item its second derivative \code{d2lpdf} (optionally);
    \item a starting partition \code{ib} of the domain of the target
    distribution such that each subinterval contains at most one
    inflection point of the transformed density;
    \item the parameter(s) \code{cT} of the transformation either for
    the entire domain or alternatively for each of the subintervals of
    the partition.
  }

  The starting partition of the domain of the target distribution into
  non-overlapping intervals has to satisfy the following conditions:
  \itemize{
    \item The partition points must be given in ascending order
    (otherwise they are sorted anyway).
    \item The first and last entry of this vector are the boundary
    points of the domain of the distribution.
    In the case when the domain of the distribution is unbounded, the
    respective points are \code{-Inf} and \code{Inf}.
    \item Within each interval of the partition, the transformed density
    possesses at most one inflection point (including all finite
    boundary points).
    \item If a boundary point is infinite, or the density vanishes at the
    boundary point, then the transformed density must be concave near
    the corresponding boundary point and in the corresponding tail,
    respectively.
    \item If the log-density \code{lpdf} has a pole or cusp at some
    point \eqn{x}, then this must be added to the starting partition
    point. Moreover, it has to be counted as inflection point.
    Moreover, in the corresponding intervals the transformed density
    must be convex.
  }

  Argument \code{d2lpdf} is optional. If \code{d2lpdf=NULL}, then
  a variant of the method is used, that determines intervals where the
  transformed density is concave or convex without means of the second
  derivative of the log-density.
  
  Parameter \code{cT} is either a single numeric, that is, the same
  transformation 
  \ifelse{html}{\out{<i>T<sub>c</sub></i>}}{\eqn{T_c}}
  is used for all subintervals of the domain,
  or it can be set independently for each of these intervals.
  In the latter case \code{length(cT)} must be equal to the number of
  intervals, that is, equal to \code{length(ib)-1}.
  For the choice of \code{cT} the following should be taken into
  consideration:
  \itemize{
    \item \code{cT=0} (the default) is most robust against numeric
    underflow or overflow.
    \item \code{cT=-0.5} has the fastest marginal generation time.
    \item One should always use \code{cT=0} or \code{cT=-0.5} 
    for intervals that contain a point where the derivative of the
    (log-) density vanishes (e.g., an extremum). For other values of
    \code{cT}, the algorithm is less accurate.
    \item For unbounded intervals \eqn{(-\inf,a]}{(-Inf,a]} or
    \eqn{[a,\inf)}{[a,Inf)}, one has to select \code{cT} such that
    \eqn{0 \ge c_T > -1}{0 >= cT > -1}.
    \item For an interval that contains a pole at one of its boundary
    points (i.e., there the density is unbounded), one has to select
    \code{cT} such that \eqn{c_T < -1}{cT < -1} and the
    transformed density is convex.
    \item If the transformed density is concave in some interval for a
    particular value of \code{cT}, then it is concave for all smaller
    values of \code{cT}.
    %% \item For values of \code{cT} less than 0, the table size increases
    %% with decreasing values of \code{cT}.
    %% \item The number of inflection points decreases with decreasing value
    %% of \code{cT}.
  }
  
  Parameter \code{rho} is a performance parameter. It defines an upper
  bound for ratio of the area below the hat function to the area below
  the squeeze function. This parameter is an upper bound of the
  rejection constant. More importantly, it provides an approximation to
  the number of (time consuming) evalutions of the log-density
  function \code{lpdf}.
  For \code{rho=1.01}, the log-density function is evaluated once for a
  sample of 300 random points. However, values of \code{rho} close to 1
  also increase the table size and thus make the setup more expensive.

  Parameter \code{max.intervals} defines the maximal number of
  subintervals and thus the maximal table size. Putting an upper bound
  on the table size prevents the algorithm from accidentally exhausting
  all of the computer memory due to invalid input.
  It is very unlikely that one has to increase the default value.
}

\value{
  Routine \code{Tinflex.setup} returns an 
  object of class \code{"Tinflex"} that stores the random variate
  generator (density, hat and squeeze functions, cumulated areas below
  hat). For details see sources of the algorithm or execute
  \code{print(gen,debug=TRUE)} with an object \code{gen} of class
  \code{"Tinflex"}.

  Routine \code{Tinflex.setup.C} is equivalent to \code{Tinflex.setup}
  but does all computations entirely in C. It returns an object of class
  \code{"TinflexC"} which is equivalent to class \code{"Tinflex"} but
  stores all data in an C structure instead of an R list.
}

\section{Warning }{
  It is very important to note that the user is responsible for the
  correctness of the supplied arguments. Since the algorithm works (in theory)
  for all distributions with piecewise twice differentiable density
  functions, it is not possible to detect improper arguments. It is thus
  recommended that the user inspect the generator object visually by
  means of the \code{plot} method (see \code{\link{plot.Tinflex}} for
  details).
}

\author{
  Josef Leydold \email{josef.leydold@wu.ac.at}, Carsten Botts and
  Wolfgang \enc{Hörmann}{Hoermann}.
}

\references{
  C. Botts, W. \enc{Hörmann}{Hoermann}, and J. Leydold (2013),
  Transformed Density Rejection with Inflection Points,
  Statistics and Computing 23(2), 251--260,
  \doi{10.1007/s11222-011-9306-4}
  (see also 
  Research Report Series / Department of Statistics and Mathematics
  Nr. 110, Department of Statistics and Mathematics, WU Wien,
  \url{https://epub.wu.ac.at/id/eprint/3158}).
  
  W. \enc{Hörmann}{Hoermann}, and J. Leydold (2022),
  A Generalized Transformed Density Rejection Algorithm,
  in preparation.
}

\seealso{
  See \code{\link{Tinflex.sample}} for drawing random samples,
  \code{\link{plot.Tinflex}} and \code{\link{print.Tinflex}} for
  printing and plotting objects of class \code{"Tinflex"}.
}

\examples{
## Example 1: Bimodal density
##   Density  f(x) = exp( -|x|^alpha + s*|x|^beta + eps*|x|^2 )
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
gen <- Tinflex.setup.C(lpdf, dlpdf, d2lpdf, ib=c(-Inf,0,Inf), rho=1.1)

## Print data about generator object.
print(gen)

## Draw a random sample
Tinflex.sample(gen, n=10)

## Inspect hat and squeeze visually in original scale
plot(gen, from=-2.5, to=2.5)
## ... and in transformed (log) scale.
plot(gen, from=-2.5, to=2.5, is.trans=TRUE)

## With Version 2.0 the setup also works without providing the
## second derivative of the log-density
gen <- Tinflex.setup.C(lpdf, dlpdf, d2lpdf=NULL, ib=c(-Inf,0,Inf), rho=1.1)
Tinflex.sample(gen, n=10)

## -------------------------------------------------------------------
## Example 2: Exponential Power Distribution
##   Density  f(x) = exp( -|x|^alpha )
##   with alpha > 0  [ >= 0.015 due to limitations of FPA ]

alpha <- 0.5

## Log-density and its derivatives.
lpdf   <- function(x) { -abs(x)^alpha }
dlpdf  <- function(x) { if (x==0) {0} else {-sign(x) * alpha*abs(x)^(alpha-1)}}
d2lpdf <- function(x) { -alpha*(alpha-1)*abs(x)^(alpha-2) }

## Parameter cT=-0.5:
##   There are two inflection points on either side of 0 and
##   a cusp at 0. Thus we need a partition point that separates
##   the inflection points from the cusp.
ib <- c(-Inf, -(1-alpha)/2, 0, (1-alpha)/2, Inf)

## Create generator object with c = -0.5.
gen <- Tinflex.setup.C(lpdf, dlpdf, d2lpdf, ib=ib, cT=-0.5, rho=1.1)

## Print data about generator object.
print(gen)

## Draw a random sample.
Tinflex.sample(gen, n=10)

## Inspect hat and squeeze visually in original scale
plot(gen, from=-4, to=4)
## ... and in transformed (log) scale.
plot(gen, from=-4, to=4, is.trans=TRUE)

## With Version 2.0 the setup also works without providing the
## second derivative of the log-density
gen <- Tinflex.setup.C(lpdf, dlpdf, d2lpdf=NULL, ib=ib, cT=-0.5, rho=1.1)
Tinflex.sample(gen, n=10)

## -------------------------------------------------------------------
## Example 3: Generalized Inverse Gaussian Distribution
##   Density  f(x) = x^(lambda-1) * exp(-omega/2 * (x+1/x))   x>= 0
##   with 0 < lambda < 1 and 0 < omega <= 0.5

la <- 0.4     ## lambda
om <- 1.e-7   ## omega

## Log-density and its derivatives.
lpdf   <- function(x) { ifelse (x==0., -Inf, ((la - 1) * log(x) - om/2*(x+1/x))) }
dlpdf  <- function(x) { if (x==0) { Inf} else {(om + 2*(la-1)*x-om*x^2)/(2*x^2)} }
d2lpdf <- function(x) { if (x==0) {-Inf} else {-(om - x + la*x)/x^3} }

## Parameter cT=0 near 0 and cT=-0.5 at tail:
ib <- c(0, (3/2*om/(1-la) + 2/9*(1-la)/om), Inf)
cT <- c(0,-0.5)

## Create generator object.
gen <- Tinflex.setup.C(lpdf, dlpdf, d2lpdf, ib=ib, cT=cT, rho=1.1)

## Print data about generator object.
print(gen)

## Draw a random sample.
Tinflex.sample(gen, n=10)

## Inspect hat and squeeze visually in original scale
plot(gen, from=0, to=5)

## With Version 2.0 the setup also works without providing the
## second derivative of the log-density
gen <- Tinflex.setup.C(lpdf, dlpdf, d2lpdf=NULL, ib=ib, cT=cT, rho=1.1)
Tinflex.sample(gen, n=10)

}

\keyword{distribution}
\keyword{datagen}
