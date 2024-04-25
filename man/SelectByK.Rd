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
          ClusterScalar = 1,
          MaxClusters = 15L,
          ReturnAllCommunities = FALSE,
          Verbose = FALSE,
          ShowPlot = FALSE,
          RetainHighest = TRUE)
}
\arguments{
  \item{Pairs}{
An object of class \code{PairSummaries}.
}
  \item{UserConfidence}{
A numeric value greater than 0 and less than 1 that represents a minimum PID centroid that users believe represents a \code{TRUE} predicted pair.
}
  \item{ClusterScalar}{
A numeric value used to scale selection of how many clusters are used in kmeans clustering. Total within-cluster sum of squares are fit to a right hyperbola, and the half-max is used to select cluster number. ``ClusterScalar'' is multiplied by the half-max to adjust cluster number selection.
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
  \item{RetainHighest}{
Logical indicating whether to retain the cluster with the highest PID centroid in the case where the PID is below the specified user confidence.
}
}
\details{
\code{SelectByK} uses a naive k-means routine to select for predicted pairs that belong to clusters whose centroids are greater than or equal to the user specified PID confidence. This means that the confidence is not a minimum, and that pairs with PIDs below the user confidence can be retained. The sum of within cluster sum of squares is used to approximate ``knee'' selection with the user supplied ``ClusterScalar'' value. By default, with a ``ClusterScalar'' value of 1 the half-max of a right-hyperbola fitted to the sum of within-cluster sum of squares is used to pick the cluster number for evaluation, ``ClusterScalar'' is multiplied by the half-max to tune cluster number selection.
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
# this function will be deprecated soon,
# please see the new ClusterByK() function.

DBPATH <- system.file("extdata",
                      "Endosymbionts_v02.sqlite",
                      package = "SynExtend")

data("Endosymbionts_LinkedFeatures", package = "SynExtend")

Pairs <- PairSummaries(SyntenyLinks = Endosymbionts_LinkedFeatures,
                       PIDs = TRUE,
                       Score = TRUE,
                       DBPATH = DBPATH,
                       Verbose = TRUE)

Pairs02 <- SelectByK(Pairs = Pairs)
}
