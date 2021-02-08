\encoding{UTF-8}

\name{Tinflex.sample}
\alias{Tinflex.sample}
\alias{Tinflex.sample.C}

\title{
  Draw Random Sample from Tinflex Generator Object 
}

\description{
  Draw a random sample from a generator object of class
  \code{"Tinflex"} or \code{"TinflexC"}.
}

\usage{
Tinflex.sample(gen, n=1)
Tinflex.sample.C(gen, n=1)
}

\arguments{
  \item{gen}{an object of class \code{"Tinflex"} or \code{"TinflexC"}.}
  \item{n}{sample size. (integer)}
}

\details{
  Routine \code{Tinflex.sample.C} allows objects of class
  \code{"TinflexC"} only and thus is a bit faster than the same call
  with routine \code{Tinflex.sample}.
}

\author{
  Josef Leydold \email{josef.leydold@wu.ac.at}, Carsten Botts and
  Wolfgang \enc{Hörmann}{Hoermann}.
}

\seealso{
  See \code{\link{Tinflex.setup}} for examples.
}

%% \examples{}

\keyword{distribution}
\keyword{datagen}
