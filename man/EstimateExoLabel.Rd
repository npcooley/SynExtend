\name{EstimateExoLabel}
\alias{EstimateExoLabel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate ExoLabel Disk Consumption
}
\description{
Estimate the total disk consumption for \code{\link{ExoLabel}}.
}
\usage{
EstimateExoLabel(num_v, avg_degree=2,
              num_edges=num_v*avg_degree,
              node_name_length=10L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{num_v}{
Approximate number of total unique nodes in the network.
}

\item{avg_degree}{
Average degree of each node in the network.
}

\item{num_edges}{
Approximate total number of edges in the network.
}

\item{node_name_length}{
Approximate average length of each node name, in characters.
}
}
\details{
This function provides a rough estimate of the total disk space required to run \code{\link{ExoLabel}} for a given input network. \code{avg_degree} and \code{num_edges} need not both be specified. The function prints out the estimated size of the original edgelist files, the estimated disk space and RAM to be consumed by \code{\link{ExoLabel}}, and the approximate ratio of disk space relative to the original file.

\code{node_name_length} specifies the average length of the node names--since the names themselves must be stored on disk, this contributes to the overall size. For relatively short node names (1-16 characters) this has a negligible impact on overall disk consumption, though it may impact the worst-case RAM consumption. Expected RAM consumption is determined by the average prefix length a random pair of vertex labels have in common, and should be closer to the minimum usage in most scenarios (see \code{\link{ExoLabel}} for more details on this).
}
\value{
Invisibly returns a vector of length six, showing the estimated RAM consumption, estimated input edgelist file size, estimated disk consumption using in-place sort (\code{use_fast_sort=FALSE}), estimated disk consumption using fast sort (\code{use_fast_sort=TRUE}), estimated final file size, and ratio of the input file size to total ExoLabel disk usage. All values denote bytes.
}

\author{
Aidan Lakshman <AHL27@pitt.edu>
}
\note{
Estimating the average node label size is challenging, and unfortunately it does have a relatively large effect on the estimated edgelist file size. This function should be used for \bold{rough} estimations of sizing, not absolute values. Errors in estimation of rough node name size will have a larger impact on edgelist file estimation than on the ExoLabel disk usage, so users can have higher confidence in estimated ExoLabel consumption.
}

\seealso{
\code{\link{ExoLabel}}
}
\examples{
# 100,000 nodes, average degree 2
EstimateExoLabel(num_v=100000, avg_degree=2)

# 10,000 nodes, 50,000 edges
EstimateExoLabel(num_v=10000, num_edges=50000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
