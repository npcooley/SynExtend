\name{ProtWeb}
\alias{ProtWeb}
\title{
ProtWeb: Predictions from ProtWeaver
}
\docType{data}
\description{
ProtWeb objects are outputted from \code{\link{predict.ProtWeaver}}. 

This class wraps the \code{\link{simMat}} object with some other diagnostic
information intended to help interpret the output of \code{\link{ProtWeaver}}
predictions..
}

\details{
\code{\link{predict.ProtWeaver}} returns a \code{ProtWeb} object, which bundles some methods
to make formatting and printing of results slightly nicer. This currently
only implements a \code{plot} function, but future functionality is in the works.
}
\value{
An object of class \code{"ProtWeb"}, which inherits from \code{"simMat"}.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\seealso{
\code{\link{predict.ProtWeaver}}

\code{\link{simMat}}

\code{\link{plot.ProtWeb}}
}
\examples{
##############
## Prediction with built-in model and data 
###############

exData <- get(data("ExampleStreptomycesData"))

# Subset isn't necessary but is faster for a working example
pw <- ProtWeaver(exData$Genes[1:10])

protweb <- predict(pw, method='Jaccard')

# print out results as an adjacency matrix
print(protweb)

# print out results as a pairwise data.frame
as.data.frame(protweb)
}