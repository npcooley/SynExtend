\name{ExampleStreptomycesData}
\alias{ExampleStreptomycesData}
\docType{data}
\title{
Example EvoWeaver Input Data from \emph{Streptomyces} Species
}
\description{
Data from \emph{Streptomyces} species to test \code{\link{EvoWeaver}} functionality.
}
\usage{data("ExampleStreptomycesData")}
\value{
  The data contain two elements, \code{Genes} and \code{Tree}. \code{Genes} is
  a list of presence/absence vectors in the input required for \code{\link{EvoWeaver}}.
  \code{Tree} is a species tree used for additional input.
}
\details{
This dataset contains a number of Clusters of Orthologous Genes (COGs) and a species
tree for use with EvoWeaver. This dataset showcases an example of using EvoWeaver with
a list of vectors. Entries in each vector are formatted correctly for use with
co-localization prediction. Each COG \code{i} contains entries of the form
\code{a_b_c}, indicating that the gene was found in genome \code{a} on chromosome
\code{b}, and was at the \code{c}'th location. The original dataset is comprised
of 301 unique genomes.
}
\seealso{
\code{\link{EvoWeaver}}
}
\examples{
exData <- get(data("ExampleStreptomycesData"))
ew <- EvoWeaver(exData$Genes)
# Subset isn't necessary but is faster for a working example
predict(ew, Subset=1:10, MySpeciesTree=exData$Tree)
}
\keyword{datasets}
