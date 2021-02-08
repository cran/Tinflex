\encoding{UTF-8}

\name{plot.Tinflex}
\alias{plot.Tinflex}
\alias{plot.TinflexC}

\title{
  Plot Tinflex Generator Objects
}

\description{
  Plotting methods for generator objects of classes
  \code{"Tinflex"} and \code{"TinflexC"}.
  The plot shows the (transformed) density, hat and squeeze.
}

\usage{
\S3method{plot}{Tinflex}(x, from, to, is.trans=FALSE, n=501, \dots)
\S3method{plot}{TinflexC}(x, from, to, is.trans=FALSE, n=501, \dots)
}

\arguments{
  \item{x}{an object of class \code{"Tinflex"} or \code{"TinflexC"}.}
  \item{from,to}{the range over which the function will be plotted. (numeric)}
  \item{is.trans}{if \code{TRUE} then the transformed density and its
    hat and squeezes are plotted. (logical)}
  \item{n}{the number of x values at which (transformed) PDF to evaluate. (integer)}
  \item{\dots}{arguments to be passed to methods, such as graphical
    parameters (see \code{\link{par}}). In particular the following
    argument may be useful:
    \describe{
      \item{\code{ylim}}{limit for the plot range: see
	\code{\link{plot.window}}. It has sensible defaults if omitted.}
    }
  }
}

\details{
  This is the \code{\link{print}} method for objects of class
  \code{"Tinflex"} or \code{"TinflexC"}.
  It plots the given density function (blue) in the
  domain (\code{from},\code{to}) as well as hat function (red) and
  squeeze (green) of the acceptance-rejection algorithm.
  If \code{is.trans} is set to \code{TRUE}, then density function, hat
  and squeeze are plotted on the transformed scale.
  Notice that the latter only gives a sensible picture if parameter
  \code{cT} is the same for all intervals.
}

\author{
  Josef Leydold \email{josef.leydold@wu.ac.at}, Carsten Botts and
  Wolfgang \enc{Hörmann}{Hoermann}.
}

\seealso{
  \code{\link{plot}}, \code{\link{plot.function}}.
  See \code{\link{Tinflex.setup}} for examples.
}

%% \examples{}

\keyword{distribution}
\keyword{datagen}
