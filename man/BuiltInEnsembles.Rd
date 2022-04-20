\name{BuiltInEnsembles}
\alias{BuiltInEnsembles}
\docType{data}
\title{
Pretrained Ensemble ProtWeaver models
}
\description{
ProtWeaver has best performance with an ensemble method combining individual 
evidence streams. Users can train their own models using:

\code{data <- predict.ProtWeaver(pw, 'Ensemble',
                                noPrediction=TRUE)}

\code{myModel <- model.fit(data)}

\code{predict.ProtWeaver(pw, 'Ensemble', pretrainedmodel=myModel)}

However, pretrained models 
are provided for ease of use. All provided models are trained on genes from \emph{Streptomyces} species.
}
\usage{data(BuiltInEnsembles)}
\format{
  The format is a list of objects of class ``glm''.
}
\details{
A model for predicting the probability of functional association for a pair of
genes or proteins in a ProtWeaver object
}
\examples{
data(BuiltInEnsembles)
}
\keyword{datasets}
