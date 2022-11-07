\name{KFDist}
\alias{KFDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kuhner-Felsenstein Distance
}
\description{
Calculate KF distance between two unrooted phylogenies.
}
\usage{
KFDist(dend1, dend2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dend1}{
  An object of class \code{dendrogram}, representing an unrooted bifurcating phylogenetic
  tree.
}
  \item{dend2}{
  An object of class \code{dendrogram}, representing an unrooted bifurcating phylogenetic 
  tree.
}
}
\details{
This function implements Kuhner-Felsenstein distance between two unrooted phylogenies. Each dendrogram is first pruned to only internal branches implying a partition in the shared leaf set; trivial partitions (where one leaf set contains 1 or 0 leaves) are ignored. The total score is calculated as the sum of squared differences between lengths of branches implying equivalent partitions. If a particular branch is unique to a given tree, it is treated as having length 0 in the other tree. The final score is normalized by the sum of squared lengths of all internal branches of both trees, resulting in a final distance that ranges from 0 to 1.
}
\value{
Returns a normalized distance, with 0 indicating identical trees and 1 indicating maximal difference.

If the trees have no leaves in common, the function will return \code{1}.
}
\references{
Robinson, D.F. and Foulds, L.R. \emph{Comparison of phylogenetic trees.} Mathematical Biosciences, 1987. \bold{53}(1–2): 131–147.

Kuhner, M. K. and Felsenstein, J. \emph{Simulation comparison of phylogeny algorithms under equal and unequal evolutionary rates.} Molecular Biology and Evolution, 1994. \bold{11}: 459–468.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}
\note{
Note that this function requires the input dendrograms to be labeled alike (ex.
leaf labeled \code{abc} in \code{dend1} represents the same species as 
leaf labeled \code{abc} in \code{dend2}). 
Labels can easily be modified using \code{\link{dendrapply}}. 
}

\examples{
# making some toy dendrograms
set.seed(123)
dm1 <- as.dist(matrix(runif(64, 0.5, 5), ncol=8))
dm2 <- as.dist(matrix(runif(64, 0.5, 5), ncol=8))
names(dm1) <- names(dm2) <- c('a', 'b', 'c', 'd', 
                              'e', 'f', 'g', 'h')
tree1 <- as.dendrogram(hclust(dm1))
tree2 <- as.dendrogram(hclust(dm2))

# get KF distance
KFDist(tree1, tree2)
}

