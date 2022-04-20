\name{predict.ProtWeaver}
\alias{predict.ProtWeaver}
\title{
Make predictions with ProtWeaver objects
}
\description{
The description will go here!.
}
\usage{
\method{predict}{ProtWeaver}(object, 
method='Ensemble', subset=NULL, verbose=TRUE, returnRawData=FALSE, ...)
}
\arguments{
  \item{object}{
A ProtWeaver object
}
  \item{method}{
Method to use for prediction. See 'Details'.
}
  \item{subset}{
Subset of data to predict on.
}
  \item{verbose}{
Should output be printed to the console?
}
\item{returnRawData}{
Internal parameter used for ensemble predictions. If TRUE, returns predictions
without formatting them into a ProtWeb object.
}
\item{...}{
Other stuff
}
}
\details{
Details will go here!
}
\value{
Returns a ProtWeb object.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\seealso{
\code{\link{ProtWeaver}}
}
\examples{
# This will be an example
}
