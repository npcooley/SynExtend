\name{GlobalSelect}
\alias{GlobalSelect}
\docType{data}
\title{
Model for identifying erroneously linked pairs
}
\description{
Though the function PairSummaries provides an argument allowing users to ask for alignments, given the time consuming nature of that process on large data, models are provided which allow for the quick and efficient identification of pairs whose PID would likely fall within a random distribution of PIDs.
}
\usage{data("GlobalSelect")}
\format{
  The format is an object of class ``glm''.
}
\details{
A model for rejecting identified pairs whose link statistics indicate a likely global PID that would fall within a random distribution in an amino acid alignment.
}
\examples{
data(GlobalSelect)
}
\keyword{datasets}
