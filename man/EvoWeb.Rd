\name{EvoWeb}
\alias{EvoWeb}
\title{
EvoWeb: Predictions from EvoWeaver
}
\docType{data}
\description{
EvoWeb objects are outputted from \code{\link{predict.EvoWeaver}}.

This class wraps the \code{\link{simMat}} object with some other diagnostic
information intended to help interpret the output of \code{\link{EvoWeaver}}
predictions..
}

\details{
\code{\link{predict.EvoWeaver}} returns a \code{EvoWeb} object, which bundles some methods
to make formatting and printing of results slightly nicer. This currently
only implements a \code{plot} function, but future functionality is in the works.
}
\format{
An object of class \code{"EvoWeb"}, which inherits from \code{"simMat"}.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\seealso{
\code{\link{predict.EvoWeaver}}

\code{\link{simMat}}

\code{\link{plot.EvoWeb}}
}
\examples{
##############
## Prediction with built-in model and data
###############

exData <- get(data("ExampleStreptomycesData"))

# Subset isn't necessary but is faster for a working example
ew <- EvoWeaver(exData$Genes[1:10])

evoweb <- predict(ew, Method='ExtantJaccard')

# print out results as an adjacency matrix
print(evoweb)

# print out results as a pairwise data.frame
as.data.frame(evoweb)
}
