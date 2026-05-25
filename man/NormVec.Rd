\name{NormVec}
\alias{NormVec}
\title{
Unit normalize a vector
}
\description{
This function is designed to work internally to functions within \code{SynExtend} so it works on relatively simple atomic vectors and has little overhead checking.
}
\usage{
NormVec(vec)
}
\arguments{
  \item{vec}{
  A numeric or integer vector.
}
}
\details{
  \code{NormVec} unit normalized a vector.
}
\value{
  A numeric vector the same length as the input.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{NucleotideOverlap}}, \code{\link{SummarizePairs}}
}
\examples{
x <- NormVec(rnorm(n = 50, mean = 2, sd = 2))
}
