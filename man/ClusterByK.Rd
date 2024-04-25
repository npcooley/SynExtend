\name{ClusterByK}
\alias{ClusterByK}
\title{
Predicted pair trimming using K-means.
}
\description{
A relatively simple k-means clustering approach to drop predicted pairs that belong to clusters with a PID centroid below a specified user threshold.
}
\usage{
ClusterByK(SynExtendObject,
           UserConfidence = list("PID" = 0.3),
           ClusterScalar = 4,
           MaxClusters = 15L,
           ColSelect = c("p1featurelength",
                         "p2featurelength",
                         "TotalMatch",
                         "Consensus",
                         "PID",
                         "Score"),
           ColNorm = "Score",
           ShowPlot = FALSE,
           Verbose = FALSE)
}
\arguments{
  \item{SynExtendObject}{
An object of class \code{PairSummaries}.
}
  \item{UserConfidence}{
A named list of length 1 where the name identifies a column of the \code{PairSummaries} object, and the value identifies a user confidence. Every k-means cluster with a center value of the column value selected greater than the confidence is retained.
}
  \item{ClusterScalar}{
A numeric value used to scale selection of how many clusters are used in kmeans clustering. A transformed total within-cluster sum of squares value is fit to a right hyperbola, and a scaled half-max value is used to select cluster number. ``ClusterScalar'' is multiplied by the half-max to adjust cluster number selection.
}
  \item{MaxClusters}{
Integer value indicating the largest number of clusters to test in a series of k-means clustering tests.
}

  \item{ColSelect}{
A character vector of column names indicating which columns to use for k-means clustering. When ``p1featurelength'', ``p2featurelength'', and ``TotalMatch'' are included together, they are morphed into a value representing the match size proportional to the longer of the two sequences.
}
  \item{ColNorm}{
A character vector of column names indicating columns the user would like to unit normalize. By default only set to ``Score''.
}
  \item{ShowPlot}{
Logical indicating whether or not to plot the CDFs for the PIDs of all k-means clusters for the determined cluster number.
}
  \item{Verbose}{
Logical indicating whether or not to display a progress bar and print the time difference upon completion.
}
}
\details{
\code{ClusterByK} uses a naive k-means routine to select for predicted pairs that belong to clusters whose centroids are greater than or equal to the user specified column-value pair. This means that the confidence is not a minimum, and that pairs with values below the user confidence can be retained. The sum of within cluster sum of squares is used to approximate ``knee'' selection with the ``ClusterScalar'' value. With a ``ClusterScalar'' value of 1 the half-max of a right-hyperbola fitted to the sum of within-cluster sum of squares is used to pick the cluster number for evaluation, ``ClusterScalar'' is multiplied by the half-max to tune cluster number selection.
\code{ClusterByK} returns the original object with an appended column and new attributes. The new column ``ClusterID'' is an integer value indicating which k-means cluster a candidate pair belongs to, while the attribute ``Retain'' is a named logical vector where the names correspond to ClusterIDs, and the logical value indicates whether the cluster center was above the user suppled column-value pair.
This function is intended to be used at the genome-to-genome comparison level, and not say, at the level of an all-vs-all comparison of many genomes. It will work well in all-vs-all cases, but it is not optimized for that scale yet.
}
\value{
An object of class \code{PairSummaries}.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{SummarizePairs}}, \code{\link{NucleotideOverlap}}, \code{\link{FindSynteny}}, \code{\link{ExpandDiagonal}}
}
\examples{
data("Endosymbionts_Pairs01", package = "SynExtend")

Pairs02 <- ClusterByK(SynExtendObject = Endosymbionts_Pairs01)
}
