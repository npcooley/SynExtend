\name{BuiltInEnsembles}
\alias{BuiltInEnsembles}
\docType{data}
\title{
Pretrained EvoWeaver Ensemble Models
}
\description{
EvoWeaver has best performance with an ensemble method combining individual
evidence streams.
This data file provides pretrained models for ease of use.
Two groups of models are provided:
1. Models trained on the KEGG MODULES dataset
2. Models trained on the CORUM dataset

These models are used internally if the user does not provide their own model,
and aren't explicitly designed to be accessed by the user.

See the examples for how to train your own ensemble model.
}
\usage{data("BuiltInEnsembles")}
\value{
  The data contain a named list of objects of class \code{\link[stats]{glm}}.
  This list currently has two entries: \code{"KEGG"} and \code{"CORUM"}
}

\examples{
## Training own ensemble method to avoid
## using built-ins

exData <- get(data("ExampleStreptomycesData"))
ew <- EvoWeaver(exData$Genes[seq_len(50L)], MySpeciesTree=exData$Tree)
datavals <- predict(ew, NoPrediction=TRUE, Verbose=interactive())

# Make sure the actual values correspond to the right pairs!
# This example just picks random numbers
# Do not do this for your own models
actual_values <- sample(c(0,1), nrow(datavals), replace=TRUE)
datavals[,'y'] <- actual_values
myModel <- glm(y~., datavals[,-c(1,2)], family='binomial')

predictionPW <- EvoWeaver(exData$Genes[51:60], MySpeciesTree=exData$Tree)
predict(predictionPW,
          PretrainedModel=myModel, Verbose=interactive())
}
\keyword{datasets}
