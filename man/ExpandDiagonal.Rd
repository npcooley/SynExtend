\name{ExpandDiagonal}
\alias{ExpandDiagonal}
\title{
Attempt to expand blocks of paired features in a \code{PairSummaries} object.
}
\description{
Attempt to expand blocks of paired features in a \code{PairSummaries} object.
}
\usage{
ExpandDiagonal(SynExtendObject,
               FeatureSeqs,
               DataBase,
               InheritConfidence = TRUE,
               GapTolerance = 100L,
               DropSingletons = FALSE,
               UserConfidence = list("PID" = 0.3),
               Verbose = FALSE)
}
\arguments{
  \item{SynExtendObject}{
An object of class \code{PairSummaries}.
}
  \item{FeatureSeqs}{
An object of class \code{FeatureSeqs}.
}
  \item{DataBase}{
A character string pointing to a SQLite database, or a connection to a \code{DECIPHER} database.
}
  \item{InheritConfidence}{
A logical indicating whether or not to inheret the user specified column-value pairs assigned to the input object.
}
  \item{GapTolerance}{
Integer value indicating the \code{diff} between feature IDs that can be tolerated to view features as part of the same block. Set by default to \code{100L}.
}
  \item{DropSingletons}{
Ignore solo pairs when planning expansion routes. Set to \code{FALSE} by default.
}
  \item{UserConfidence}{
A named list of length 1 where the name identifies a column of the \code{PairSummaries} object, and the value identifies a user confidence. Every k-means cluster with a center value of the column value selected greater than the confidence is retained.
}
  \item{Verbose}{
Logical indicating whether or not to display a progress bar and print the time difference upon completion.
}
}
\details{
\code{ExpandDiagonal} uses a naive expansion algorithm to attempt to fill in gaps in blocks of paired features and to attempt to expand blocks of paired features.
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
DBPATH <- system.file("extdata",
                      "Endosymbionts_v02.sqlite",
                      package = "SynExtend")

data("Endosymbionts_LinkedFeatures", package = "SynExtend")
Endosymbiont_Seqs <- PrepareSeqs(SynExtendObject = Endosymbionts_LinkedFeatures,
                                 DataBase = DBPATH,
                                 Verbose = TRUE)
                      
data("Endosymbionts_Pairs02", package = "SynExtend")
Pairs03 <- ExpandDiagonal(SynExtendObject = Endosymbionts_Pairs02,
                          DataBase = DBPATH,
                          FeatureSeqs = Endosymbiont_Seqs,
                          Verbose = TRUE)
}
