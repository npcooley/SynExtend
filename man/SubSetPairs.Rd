\name{SubSetPairs}
\alias{SubSetPairs}
\title{
Subset a ``PairSummaries'' object.
}
\description{
For a given object of class ``PairSummaries'', pairs based on either competing predictions, user thresholds on prediction statistics, or both.
}
\usage{
SubSetPairs(CurrentPairs,
            UserThresholds,
            RejectCompetitors = TRUE,
            RejectionCriteria = "PID",
            WinnersOnly = TRUE,
            Verbose = FALSE)
}
\arguments{
  \item{CurrentPairs}{
An object of class ``PairSummaries''. Can also take in a generic ``data.frame'', as long as the feature naming scheme is the same as that followed by all \code{SynExtend} functions.
}
  \item{UserThresholds}{
A named vector where values indicate a threshold for statistics to be above, and names designate which statistic to threshold on.
}
  \item{RejectCompetitors}{
A logical that defaults to ``TRUE''. Allowing users to choose to remove competing predictions. When set to ``FALSE'', no competitor rejection is performed. When ``TRUE'' all competing pairs with the exception of the best pair as determined by ``RejectionCriteria'' are rejected. Can additionally be set to a numeric or integer, in which case only competing predictions below that value are dropped.
}
  \item{RejectionCriteria}{
A character indicating which column value competitor rejection should reference. Defaults to ``PID''.
}
  \item{WinnersOnly}{
A logical indicating whether or not to return just the pairs that are selected. Defaults to ``TRUE'' to return  a subset object of class ``PairSummaries''. When ``FALSE'', function returns a list of two ``PairSummaries'' objects, one of the selected pairs, and the second of the rejected pairs.
}
  \item{Verbose}{
Logical indicating whether or not to display a progress bar and print the time difference upon completion.
}
}
\details{
\code{SubSetPairs} uses a naive competitor rejection algorithm to remove predicted pairs when nodes are predicted to be paired to multiple nodes within the same index.
}
\value{
An object of class ``PairSummaries'', or a list of two ``PairSummaries'' objects.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{PairSummaries}} \code{\link{NucleotideOverlap}}
}
\examples{
\dontrun{
DBPATH <- system.file("extdata",
                      "VignetteSeqs.sqlite",
                      package = "SynExtend")
Syn <- FindSynteny(dbFile = DBPATH)
GeneCalls <- vector(mode = "list",
                    length = ncol(Syn))

GeneCalls[[1L]] <- gffToDataFrame(GFF = system.file("extdata",
                                                    "GCA_006740685.1_ASM674068v1_genomic.gff.gz",
                                                    package = "SynExtend"),
                                  Verbose = TRUE)
GeneCalls[[2L]] <- gffToDataFrame(GFF = system.file("extdata",
                                                    "GCA_000956175.1_ASM95617v1_genomic.gff.gz",
                                                    package = "SynExtend"),
                                  Verbose = TRUE)
GeneCalls[[3L]] <- gffToDataFrame(GFF = system.file("extdata",
                                                    "GCA_000875775.1_ASM87577v1_genomic.gff.gz",
                                                    package = "SynExtend"),
                                  Verbose = TRUE)
names(GeneCalls) <- seq(length(GeneCalls))
Links <- NucleotideOverlap(SyntenyObject = Syn,
                           GeneCalls = GeneCalls,
                           LimitIndex = FALSE,
                           Verbose = TRUE)
PredictedPairs <- PairSummaries(SyntenyLinks = Links,
                                DBPATH = DBPATH,
                                PIDs = FALSE,
                                AcceptContigNames = TRUE,
                                Verbose = TRUE)
# remove competitors under default conditions
Pairs2 <- SubSetPairs(CurrentPairs = PredictedPairs,
                      Verbose = TRUE)
THRESH <- c(0.5, 21)
names(THRESH) <- c("Consensus", "ExactMatch")
# remove pairs only based on user defined thresholds
Pairs3 <- SubSetPairs(CurrentPairs = PredictedPairs,
                      UserThresholds = THRESH,
                      RejectCompetitors = FALSE,
                      Verbose = TRUE)
}
}
