\name{getProtWebData.ProtWeb}
\alias{getProtWebData.ProtWeb}
\alias{getProtWebData}
\title{
Extract information from a ProtWeb object
}
\description{
ProtWeb objects are outputted from \code{\link{predict.ProtWeaver}}. 

This function extracts the underlying data from the object.
}
\usage{
\method{getProtWebData}{ProtWeb}(x, asDf=FALSE, ...)
}
\arguments{
  \item{x}{
    A ProtWeb object
  }
  \item{asDf}{
    Should data be printed as a pairwise entry? If \code{TRUE}, returns a matrix 
    with three columns, where the first two columns define the pair of genes/proteins 
    and the third column defines the prediction/score. If \code{FALSE}, returns an 
    adjacency matrix encoding the same information (but may be sparse depending on
    how many predictions were made in the original \code{predict} call.)
  }
  \item{...}{
    Additional parameters for consistency with generic.
  }
}

\details{
\code{\link{predict.ProtWeaver}} returns a \code{ProtWeb} object, which bundles some methods
to make formatting and printing of results slightly nicer. This method extracts data form the \code{ProtWeb} object.

If \code{asDf=TRUE}, the return data is a \code{3xN} \code{\link[base]{data.frame}}, with columns 3 showing
the prediction for the pair of genes/proteins specified in columns 1 and 2. This format 
is when predictions are made on small number of pairs (meaning the resulting adjacency 
matrix is sparse).

If \code{asDf=FALSE}, the return data is a \code{NxN} adjacency matrix, with entry
\code{i,j} containing the prediction for genes \code{i} and \code{j}. 
}
\value{
Either a \code{\link[base]{data.frame}} or a matrix.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\seealso{
\code{\link{predict.ProtWeaver}}
}
\examples{
##############
## Prediction with built-in model and data 
###############

exData <- get(data("ExampleStreptomycesData"))
pw <- ProtWeaver(exData$Genes)

# Subset isn't necessary but is faster for a working example
protweb <- predict(pw, method='Jaccard', subset=1:50)

# print out results as an adjacency matrix
getProtWebData(protweb)

# print out results as a pairwise data.frame
getProtWebData(protweb, asDf=TRUE)
}