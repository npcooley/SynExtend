\name{PhyloDistance.JRFDist}
\alias{JRFDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Jaccard-Robinson-Foulds Distance
}
\description{
Calculate JRF distance between two unrooted phylogenies.
}
\usage{
# Internally called using PhyloDistance with Method arg
PhyloDistance(dend1, dend2, Method="JRF", RawScore=FALSE, ExpVal=2)
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
  If \code{TRUE}, returns the distance score and the total number of splits for each tree. If
  \code{FALSE}, returns the distance between the two trees.
}
\item{ExpVal}{
  Numeric specifying the value of \code{k} to use in calculation. Larger values penalize partial matches more; see Details for more information.
}
}
\details{
This function implements the Jaccard-Robinson-Foulds distance between two unrooted phylogenies. Each dendrogram is first pruned to only internal branches implying a partition in the shared leaf set; trivial partitions (where one leaf set contains 1 or 0 leaves) are ignored. 

The total score is calculated by pairing branches and scoring their similarity. For a set of two branches \eqn{A, B} that partition the leaves into \eqn{(A_1, A_2)} and \eqn{(B_1, B_2)} (resp.), the distance between the branches is calculated as:

\eqn{2 - 2\left(\frac{|X \cap Y|}{| X\cup Y|}\right)^k}

where \eqn{X \in (A_1, A_2),\; Y \in (B_1, B_2)} are chosen to maximize the score of the pairing, and \eqn{k} the value of \code{ExpVal}. The sum of these scores for all branches produces the overall distance between the two trees, which is then normalized by the number of branches in each tree.

There are a few special cases to this distance. If \code{ExpVal=1}, the distance is equivalent to the metric introduced in Nye et al. (2006). As \code{ExpVal} approaches infinity, the value becomes close to the (non-Generalized) Robinson Foulds Distance.
}
\value{
Returns a normalized distance, with 0 indicating identical trees and 1 indicating maximal difference.

If \code{RawScore=TRUE}, returns a named length 3 vector with the first entry the summed distance score over the branch pairings, and the subsequent entries the number of partitions for each tree.

If the trees have no leaves in common, the function will return \code{1} if 
\code{RawScore=FALSE}, and \code{c(0, NA, NA)} if \code{TRUE}.
}
\references{
Nye, T. M. W., Liò, P., & Gilks, W. R. \emph{A novel algorithm and web-based tool for comparing two alternative phylogenetic trees.} Bioinformatics, 2006. \bold{22}(1): 117–119.

Böcker, S., Canzar, S., & Klau, G. W.. \emph{The generalized Robinson-Foulds metric}. Algorithms in Bioinformatics, 2013. \bold{8126}: 156–169.
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

tree1 <- as.dendrogram(hclust(dm1))
tree2 <- as.dendrogram(hclust(dm2))

# Nye Metric
PhyloDistance(tree1, tree2, Method="JRF", ExpVal=1)

# Jaccard-RobinsonFoulds
JRFDist(tree1, tree2, Method="JRF", ExpVal=2)

# Good approximation to RF Dist (note RFDist is much faster for this)
JRFDist(tree1, tree2, Method="JRF", ExpVal=1000)
RFDist(tree1, tree2, Method="RF")
}

