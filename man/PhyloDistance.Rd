\name{PhyloDistance}
\alias{PhyloDistance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Distance between Unrooted Phylogenies
}
\description{
Calculates distance between two unrooted phylogenies using a variety of metrics.
}
\usage{
PhyloDistance(dend1, dend2,
              Method=c("CI", "RF", "KF", "JRF"), 
              RawScore=FALSE, JRFExp=2)
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
\item{Method}{
  Method to use for calculating tree distances. The following values are supported: \code{"CI", "RF", "KF", "JRF"}. See Details for more information.
}
\item{RawScore}{
  If \code{FALSE}, returns distance between the two trees. If \code{TRUE}, returns the component values used to calculate the distance. This may be preferred for methods like \code{GRF}. See the pages specific to each algorithm for more information on what values are reported.
}
\item{JRFExp}{
\code{k}-value used in calculation of JRF Distance. Unused if \code{Method} is not \code{"JRF"}.
}
}
\details{
This function implements a variety of tree distances, specified by the value of \code{Method}. The following values are supported, along with links to documentation pages for each function:

\itemize{
  \item{\code{"RF"}: \link[=RFDist]{Robinson-Foulds Distance}}
  \item{\code{"CI"}: \link[=CIDist]{Clustering Information Distance}}
  \item{\code{"JRF"}: \link[=JRFDist]{Jaccard-Robinson-Foulds Distance}, equivalent to the Nye Distance Metric when \code{JRFVal=1}}
  \item{\code{"KF"}: \link[=KFDist]{Kuhner-Felsenstein Distance}}
}

Information on each of these algorithms, how scores are calculated, and references to literature can be found at the above links. Method \code{"CI"} is selected by default due to recent work showing this method as the most robust tree distance metric under general conditions.
}
\value{
Returns a normalized distance, with 0 indicating identical trees and 1 indicating 
maximal difference. If the trees have no leaves in common, the function will return \code{1}.

If \code{RawScore=TRUE}, returns a vector of the components used to calculate the distance. This is typically a length 3 vector, but specific details can be found on the description for each algorithm.
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

\seealso{
\link[=RFDist]{Robinson-Foulds Distance}

\link[=CIDist]{Clustering Information Distance}

\link[=JRFDist]{Jaccard-Robinson-Foulds Distance}

\link[=KFDist]{Kuhner-Felsenstein Distance}
}
\examples{
# making some toy dendrograms
set.seed(123)
dm1 <- as.dist(matrix(runif(64, 0.5, 5), ncol=8))
dm2 <- as.dist(matrix(runif(64, 0.5, 5), ncol=8))

tree1 <- as.dendrogram(hclust(dm1))
tree2 <- as.dendrogram(hclust(dm2))

# Robinson-Foulds Distance
PhyloDistance(tree1, tree2, Method="RF")

# Clustering Information Distance
PhyloDistance(tree1, tree2, Method="CI")

# Kuhner-Felsenstein Distance
PhyloDistance(tree1, tree2, Method="KF")

# Nye Distance Metric
PhyloDistance(tree1, tree2, Method="JRF", JRFExp=1)

# Jaccard-Robinson-Foulds Distance
PhyloDistance(tree1, tree2, Method="JRF", JRFExp=2)
}

