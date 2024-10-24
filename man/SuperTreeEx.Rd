\name{SuperTreeEx}
\alias{SuperTreeEx}
\docType{data}
\title{
Example Dendrograms
}
\description{
A set of four dendrograms for use in
\code{\link{SuperTree}} examples.
}
\usage{data("SuperTreeEx")}
\format{
  A list with four elements, where each is a object of type \code{\link{dendrogram}}
  corresponding to a gene tree constructed from a set of 301 \emph{Streptomyces}
  genomes. Each leaf node is labeled in the form \code{A_B_C}, where \code{A} is
  a number identifying the genome, \code{B} is a number identifying the contig,
  and \code{C} is a number identifying the gene. Altogether, each label uniquely
  identifies a gene.
}
\examples{
data(SuperTreeEx, package="SynExtend")
}
\keyword{datasets}
