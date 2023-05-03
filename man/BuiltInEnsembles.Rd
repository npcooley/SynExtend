\name{BuiltInEnsembles}
\alias{BuiltInEnsembles}
\docType{data}
\title{
Pretrained ProtWeaver Ensemble Models
}
\description{
ProtWeaver has best performance with an ensemble method combining individual 
evidence streams. 
This data file provides pretrained models for ease of use. 
These models are trained on genes from \emph{Streptomyces} species. 

These models are used internally if the user does not provide their own model,
and aren't explicitly designed to be accessed by the user.

See the examples for how to train your own ensemble model.
}
\usage{data("BuiltInEnsembles")}
\value{
  The data contain a list of objects of class \code{\link[stats]{glm}}.
}

\examples{
## Training own ensemble method to avoid 
## using built-ins

exData <- get(data("ExampleStreptomycesData"))
pw <- ProtWeaver(exData$Genes[1:50])
datavals <- predict(pw, NoPrediction=TRUE)
                  
# Make sure the actual values correspond to the right pairs! 
# This example just picks random numbers
# Do not do this for your own models
actual_values <- sample(c(0,1), nrow(datavals), replace=TRUE)
datavals[,'y'] <- actual_values                  
myModel <- glm(y~., datavals[,-c(1,2)], family='binomial')

predictionPW <- ProtWeaver(exData$Genes[51:60])
predict(predictionPW, 
          PretrainedModel=myModel)
}
\keyword{datasets}
