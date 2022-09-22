\name{SelectByK}
\alias{SelectByK}
\title{
Predicted pair trimming using K-means.
}
\description{
A relatively simple k-means clustering approach to drop predicted pairs that belong to clusters with a PID centroid below a specified user threshold.
}
\usage{
SelectByK(Pairs,
          UserConfidence = 0.5,
          ClusterSelect = 0.4,
          MaxClusters = 15L,
          ReturnAllCommunities = FALSE,
          Verbose = FALSE,
          ShowPlot = FALSE)
}
\arguments{
  \item{Pairs}{
An object of class \code{PairSummaries}.
}
  \item{UserConfidence}{
A numeric value greater than 0 and less than 1 that represents a minimum PID centroid that users believe represents a \code{TRUE} predicted pair.
}
  \item{ClusterSelect}{
A numeric value greater than 0 and less than 1 that represents a user threshold to select the ``knee'' of the tested k-means clusters for cluster number selection. The value used selects for the least number of clusters that produce a within-cluster sum of squares less than or equal to the quantile value of ``ClusterSelect'' in the vector of sum of squares values for all tested clusterings.
}
  \item{MaxClusters}{
Integer value indicating the largest number of clusters to test in a series of k-means clustering tests.
}

  \item{ReturnAllCommunities}{
A logical value, if ``TRUE'', function returns of a list where the second position is a list of ``PairSummaries'' tables for each k-means cluster. By default is ``FALSE'', returning only a ``PairSummaries'' object of the retained predicted pairs.
}
  \item{ShowPlot}{
Logical indicating whether or not to plot the CDFs for the PIDs of all k-means clusters for the determined cluster number.
}
  \item{Verbose}{
Logical indicating whether or not to display a progress bar and print the time difference upon completion.
}
}
\details{
\code{SelectByK} uses a naive k-means routine to select for predicted pairs that below to clusters whose centroids are greater than or equal to the user specified PID confidence. This means that the confidence is not a minimum, and that pairs with PIDs below the user confidence can be retained. The sum of within cluster sum of squares is used to approximate ``knee'' selection with the user supplied ``ClusterSelect'' value. This value selects the smallest number of k-means clusters with total within cluster sum of squares less than the the quantile value of ``ClusterSelect'' for the vector of total sum of squares for all clustering attempts up to the ``MaxClusters'' value.
This function is intended to be used at the genome-to-genome comparison level, and not say, at the level of an all-vs-all comparison of many genomes.
}
\value{
An object of class \code{PairSummaries}.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{PairSummaries}}, \code{\link{NucleotideOverlap}}, \code{link{SubSetPairs}}, \code{\link{FindSynteny}}
}
\examples{
data("Endosymbionts_Pairs01", package = "SynExtend")
Pairs02 <- SelectByK(Pairs = Endosymbionts_Pairs01)
}
