\name{PhyloDistance-CIDist}
\alias{PhyloDistance-CI}
\alias{CIDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clustering Information Distance
}
\description{
Calculate distance between two unrooted phylogenies using mutual clustering
information of branch partitions.
}

\details{
This function is called as part of \code{\link{PhyloDistance}} and calculates tree distance using the clustering information
approach first described in Smith (2020). This function iteratively pairs
internal tree branches of a phylogeny based on their similarity, then scores 
overall similarity as the sum of these measures. The similarity score is then
converted to a distance by normalizing by the average entropy of the two trees.
This metric has been demonstrated to outperform numerous other metrics in capabilities; 
see the original publication cited in References for more information. 

Users may wish to use the actual similarity values rather than a distance metric;
the option to specify \code{RawScore=TRUE} is provided for this case. Distance
is calculated as \eqn{\frac{M - S}{M}}, where \eqn{M=\frac{1}{2}(H_1 + H_2)}, \eqn{H_i} is the entropy of the \eqn{i}'th tree, and \eqn{S} is the similarity score between them. As shown in
the original publication, this satisfies the necessary requirements to be considered
a distance metric. Setting \code{RawScore=TRUE} will instead return a vector with
\eqn{(S, H_1, H_2, p)}, where \eqn{p} is an approximation for the two sided p-value of the result based on random simulations from Smith (2020).
}
\value{
Returns a normalized distance, with 0 indicating identical trees and 1 indicating 
maximal difference. Note that branch lengths are not considered, so two trees
with different branch lengths may return a distance of 0.

If \code{RawScore=TRUE}, returns a named length 4 vector with the first entry the similarity
score, subsequent entries the entropy values for each tree, and the last entry the approximate p-value for the result based on simulations.

If the trees have no leaves in common, the function will return \code{1} if 
\code{RawScore=FALSE}, and \code{c(0, NA, NA, NA)} if \code{TRUE}.
}
\references{
Smith, Martin R. \emph{Information theoretic generalized Robinson–Foulds metrics for comparing phylogenetic trees.} Bioinformatics, 2020. \bold{36}(20):5007-5013.
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

# get RF distance
PhyloDistance(tree1, tree2, Method="CI")

# get similarity score with individual entropies
PhyloDistance(tree1, tree2, Method="CI", RawScore=TRUE)
}

