\name{CIDist_NullDist}
\alias{CIDist_NullDist}
\docType{data}
\title{
Simulated Null Distributions for CI Distance
}
\description{
Simulated values of \link[=PhyloDistance-CI]{Clustering Information Distance} for random trees with 4 to 200 shared leaves.
}
\usage{data("CIDist_NullDist")}
\format{
  A matrix \code{CI_DISTANCE_INTERNAL} with 197 columns and 13 rows.
}
\details{
Each column of the matrix corresponds to the distribution of distances between random trees with the given number of leaves. This begins at \code{CI_DISTANCE_INTERNAL[,1]} corresponding to 4 leaves, and ends at \code{CI_DISTANCE_INTERNAL[,197]} corresponding to 200 leaves. Distances begin at 4 leaves since there is only one unrooted tree with 1, 2, or 3 leaves (so the distance between any given tree with less than 4 leaves is always 0).

Each row of the matrix corresponds to statistics for the given simulation set. The first row gives the minimum value, the next 9 give quantiles in \code{c(1\%, 5\%, 10\%, 25\%, 50\%, 75\%, 90\%, 95\%, 99\%)}, and the last three rows give the \code{max}, \code{mean}, and \code{sd} (resp.).
}
\source{
Datafiles obtained from the \href{https://ms609.github.io/TreeDistData/index.html}{TreeDistData} package, published as part of Smith (2020).
}
\references{
Smith, Martin R. \emph{Information theoretic generalized Robinson–Foulds metrics for comparing phylogenetic trees.} Bioinformatics, 2020. \bold{36}(20):5007-5013.
}
\examples{
data(CIDist_NullDist)
}
\keyword{datasets}