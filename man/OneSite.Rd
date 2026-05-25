\name{OneSite}
\alias{OneSite}
\title{
Calculate a site on a right hyperbola.
}
\description{
This function is designed to work internally to functions within \code{SynExtend} so it works on relatively simple atomic vectors and has little overhead checking.
}
\usage{
OneSite(X,
        Bmax,
        Kd)
}
\arguments{
  \item{X}{
  Numeric; an x coordinate value.
}
  \item{Bmax}{
  Numeric; an asymptotic value.
}
  \item{Kd}{
  Numeric; the half-max of the right hyperbola.
}
}
\details{
  \code{OneSite} calculates the Y-value for a given X-value on a right hyperbola.
}
\value{
  A numeric of length 1.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{NucleotideOverlap}}, \code{\link{SummarizePairs}}
}
\examples{
x <- OneSite(X = 3,
             Bmax = 10,
             Kd = 3)
             
# plot(x = 1:10, y = vapply(X = 1:10, FUN = function(x) {OneSite(X = x, Bmax = 5, Kd = 2)}, FUN.VALUE = vector(mode = "numeric", length = 1)))
}
