\name{Generic}
\alias{Generic}
\docType{data}
\title{
Model for predicting PID based on k-mer statistics
}
\description{
Though the function PairSummaries provides an argument allowing users to ask for alignments, given the time consuming nature of that process on large data, models are provided for predicting PIDs of pairs based on k-mer statistics without performing alignments.
}
\usage{data("Generic")}
\value{
  The format is an object of class ``glm''.
}
\details{
A model for predicting the PID of a pair of sequences based on the k-mers that were used to link the pair.
}
\examples{
data(Generic)
}
\keyword{datasets}
