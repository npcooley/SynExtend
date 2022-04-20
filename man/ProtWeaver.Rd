\name{ProtWeaver}
\alias{ProtWeaver}
\title{
Build functional association networks from phylogenetic data
}
\description{
The description will go here!.
}
\usage{
ProtWeaver(listOfData)
}
\arguments{
  \item{listOfData}{
A list of gene data, where each entry corresponds to information on a particular
gene. List must contain either dendrograms or vectors, and cannot contain a mixture. 
If list is composed of dendrograms, each dendrogram is a gene tree for the corresponding entry. If list is composed of vectors, vectors should be numeric or character vectors denoting the genomes containing that gene. 
}
}
\details{
Details will go here!
}
\value{
Returns a ProtWeaver object.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\seealso{
\code{\link{predict.ProtWeaver}}
}
\examples{
# This will be an example
}
