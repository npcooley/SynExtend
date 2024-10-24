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
               DataBase01,
               InheritConfidence = FALSE,
               GapTolerance = 100L,
               DropSingletons = FALSE,
               UserConfidence = list("PID" = 0.3),
               Processors = 1,
               Verbose = FALSE)
}
\arguments{
  \item{SynExtendObject}{
An object of class \code{PairSummaries}.
}
  \item{DataBase01}{
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
A named list of length 1 where the name identifies a column of the \code{PairSummaries} object, and the value identifies a user confidence. To be retained, a pair evaluated for expansion must be above all user specified confidences.
}
  \item{Processors}{
An integer value indicating how many processors to supply to \code{\link{AlignPairs}}. Supplying \code{NULL} will cause detection and use of all available cores.
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
library(RSQLite)
DBPATH <- system.file("extdata",
                      "Endosymbionts_v02.sqlite",
                      package = "SynExtend")
tmp <- tempfile()
system(command = paste("cp",
                       DBPATH,
                       tmp))
DBCONN <- dbConnect(SQLite(), tmp)
                      
data("Endosymbionts_LinkedFeatures", package = "SynExtend")
PrepareSeqs(SynExtendObject = Endosymbionts_LinkedFeatures,
            DataBase01 = DBCONN,
            Verbose = TRUE)
SummarizedPairs <- SummarizePairs(SynExtendObject = Endosymbionts_LinkedFeatures,
                                  DataBase01 = DBCONN,
                                  Verbose = TRUE)
ExpandedPairs <- ExpandDiagonal(SynExtendObject = SummarizedPairs,
                                DataBase01 = DBCONN,
                                Verbose = TRUE)
dbDisconnect(DBCONN)
unlink(tmp)

}
