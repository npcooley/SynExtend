\name{RFDist}
\alias{RFDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robinson-Foulds Distance
}
\description{
Calculate RF distance between two unrooted phylogenie..
}
\usage{
RFDist(dend1, dend2, RawScore=FALSE)
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
\item{RawScore}{
  If \code{TRUE}, returns the number of unique splits and the total number of splits for each tree. If
  \code{FALSE}, returns the distance between the two trees.
}
}
\details{
This function implements Robinson-Foulds distance between two unrooted phylogenies. Each dendrogram is first pruned to only internal branches implying a partition in the shared leaf set; trivial partitions (where one leaf set contains 1 or 0 leaves) are ignored. The total score is calculated as the number of unique partitions divided by the total number of partitions in both trees. Setting \code{RawScore=TRUE} will instead return a vector with
\eqn{(P_{shared}, P_1, P_2)}, corresponding to the shared partitions and partitions in the first and second trees (respectively).

This algorithm incorporates some optimizations from Pattengale et al. (2007) to improve computation time of the original fast RF algorithm detailed in Day (1985).
}
\value{
Returns a normalized distance, with 0 indicating identical trees and 1 indicating 
maximal difference. Note that branch lengths are not considered, so two trees
with different branch lengths may return a distance of 0.

If \code{RawScore=TRUE}, returns a named length 3 vector with the first entry the number of unique partitions, and the subsequent entries the number of partitions for each tree.

If the trees have no leaves in common, the function will return \code{1} if 
\code{RawScore=FALSE}, and \code{c(0, NA, NA)} if \code{TRUE}.
}
\references{
Robinson, D.F. and Foulds, L.R. \emph{Comparison of phylogenetic trees.} Mathematical Biosciences, 1987. \bold{53}(1–2): 131–147.

Day, William H.E. \emph{Optimal algorithms for comparing trees with labeled leaves.} Journal of classification, 1985. \bold{2}(1): 7-28.

Pattengale, N.D., Gottlieb, E.J., and Moret, B.M. \emph{Efficiently computing the Robinson-Foulds metric.} Journal of computational biology, 2007. \bold{14}(6): 724-735.
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

# get RF distance
RFDist(tree1, tree2)

# get number of unique splits per tree
RFDist(tree1, tree2, RawScore=TRUE)
}

