\name{ExampleStreptomycesData}
\alias{ExampleStreptomycesData}
\docType{data}
\title{
Example ProtWeaver Input Data from \emph{Streptomyces} Species
}
\description{
Data from {Streptomyces} species to test \code{\link{ProtWeaver}} functionality.
}
\usage{data("ExampleStreptomycesData")}
\format{
  The data contain two elements, \code{Genes} and \code{Tree}. \code{Genes} is 
  a list of presence/absence vectors in the input required for \code{\link{ProtWeaver}}.
  \code{Tree} is a species tree used for additional input.
}
\details{
This dataset contains a number of Clusters of Orthologous Genes (COGs) and a species
tree for use with ProtWeaver. This dataset showcases an example of using ProtWeaver with
a list of vectors. Entries in each vector are formatted correctly for use with
co-localization prediction. Each COG \code{i} contains entries of the form 
\code{a_b_c}, indicating that the gene was found in genome \code{a} on chromosome
\code{b}, and was at the \code{c}'th location. The original dataset is comprised 
of 301 unique genomes.
}
\seealso{
\code{\link{ProtWeaver}}
}
\examples{
exData <- get(data("ExampleStreptomycesData"))
pw <- ProtWeaver(exData$Genes)
# Subset isn't necessary but is faster for a working example
predict(pw, subset=1:10, mySpeciesTree=exData$Tree)
}
\keyword{datasets}
