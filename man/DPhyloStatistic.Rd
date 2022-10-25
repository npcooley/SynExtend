\name{DPhyloStatistic}
\alias{DPhyloStatistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
D-Statistic for Binary States on a Phylogeny
}
\description{
Calculates if a presence/absence pattern is random, Brownian, or neither with respect to a given phylogeny.
}
\usage{
DPhyloStatistic(dend, PAProfile, NumIter = 1000L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dend}{
  An object of class \code{\link{dendrogram}}
}
  \item{PAProfile}{
  A vector representing presence/absence of binary traits. See Details for more information.
}
  \item{NumIter}{
  Number of iterations to simulate for random permutation analysis.
}
}
\details{
This function implements the D-Statistic for binary traits on a phylogeny, as introduced in Fritz and Purvis (2009). The statstic is the following ratio:
\eqn{\frac{D_{obs} - D_b}{D_r - D_b}}.

Here \eqn{D_{obs}} is the D value for the input data, \eqn{D_b} is the value under simulated Brownian evolution, and \eqn{D_r} is the value under random permutation of the input data. The D value measures the sum of sister clade differences in a phylogeny weighted by branch lengths. A score close to 1 indicates phylogenetically random distribution, and a score close to 0 indicates the trait likely evolved under Brownian motion. Scores can fall outside this range; these scores are only intended as benchmark points on the scale. See the original paper cited in References for more information.

The input \code{PAProfile} supports a number of formatting options:
\itemize{
  \item{Character vector, where each element is a label of the dendrogram. Presence in the character vector indicates presence of the trait in the corresponding label.}
  \item{Integer vector of length equivalent to the number of leaves, comprised of 0s and 1s. 0 indicates absence in the corresponding leaf, and 1 indicates presence.}
  \item{Logical vector of length equivalent to number of leaves. \code{FALSE} indicates absence in the corresponding leaf, and \code{TRUE} indicates presence.}
  }
  
See Examples for a demonstration of each case.
}
\value{
Returns a numerical value. Values close to 0 indicate random distribution, and values close to 1 indicate a Brownian distribution.
}
\references{
Fritz S.A. and Purvis A. \emph{Selectivity in Mammalian Extinction Risk and Threat Types: a New Measure of Phylogenetic Signal Strength in Binary Traits.} Conservation Biology, 2010. \bold{24}(4):1042-1051.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\examples{
##########################################################
### Replicating results from Table 1 in original paper ###
##########################################################

# creates a dendrogram with 16 leaves and branch lengths all 1
distMat <- suppressWarnings(matrix(1:17, nrow=16, ncol=16))
testDend <- as.dendrogram(hclust(as.dist(distMat)))
testDend <- dendrapply(testDend, \(x){
                      attr(x, 'height') <- attr(x, 'height') / 2
                      return(x)
                    })
attr(testDend[[1]], 'height') <- attr(testDend[[2]], 'height') <- 3
attr(testDend, 'height') <- 4
plot(testDend)

set.seed(123)

# extremely clumped (should be close to -2.4)
DPhyloStatistic(testDend, as.character(1:8))

# clumped Brownian (should be close to 0)
DPhyloStatistic(testDend, as.character(c(1,2,5,6,10,12,13,14)))

# random (should be close to 1.0)
DPhyloStatistic(testDend, as.character(c(1,4:6,10,13,14,16)))

# overdispersed (should be close to 1.9)
DPhyloStatistic(testDend, as.character(seq(2,16,by=2)))

###########################################
### Different ways to create PAProfiles ###
###########################################

allLabs <- as.character(labels(testDend))

# All these ways create a PAProfile with 
# presence in members 1:4 
# and absence in members 5:16

# numeric vector:
c(rep(1,4), rep(0, length(allLabs)-4))

# logical vector:
c(rep(TRUE,4), rep(FALSE, length(allLabs)-4))

# character vector:
allLabs[1:4]
}