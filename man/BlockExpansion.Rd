\name{BlockExpansion}
\alias{BlockExpansion}
\title{
Attempt to expand blocks of paired features in a \code{PairSummaries} object.
}
\description{
Attempt to expand blocks of paired features in a \code{PairSummaries} object.
}
\usage{
BlockExpansion(PairsObject,
               GapTolerance = 2L,
               DropSingletons = FALSE,
               Floor = 0.5,
               Verbose = FALSE,
               NewPairsOnly = TRUE,
               DBPATH)
}
\arguments{
  \item{PairsObject}{
An object of class \code{PairSummaries}.
}
  \item{GapTolerance}{
Integer value indicating the \code{diff} between feature IDs that can be tolerated to view features as part of the same block. Set by default to \code{2L}, implying that a single feature missing in a run of pairs will not cause the block to be split. Setting to \code{3L} would imply that a \code{diff} of 3 between features, or a gap of 2 features, can be viewed as those features being part of the same block.
}
  \item{DropSingletons}{
Ignore solo pairs when planning expansion routes. Set to \code{FALSE} by default.
}
  \item{Floor}{
Lower PID limit for keeping a pair that was evaluated during expansion.
}
  \item{Verbose}{
Logical indicating whether or not to display a progress bar and print the time difference upon completion.
}
  \item{NewPairsOnly}{
Logical indicating whether or not to return only the pairs that were kept from all expansion attempts, or to return a \code{PairSummaries} object with the new pairs folded in.
}
  \item{DBPATH}{
A file or connection pointing to the \code{DECIPHER} database supplied to \code{FindSynteny} for the original map construction.
}
}
\details{
\code{BlockExpansion} uses a naive expansion algorithm to attempt to fill in gaps in blocks of paired features and to attempt to expand blocks of paired features.
}
\value{
An object of class \code{PairSummaries}.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{PairSummaries}} \code{\link{NucleotideOverlap}} \code{link{SubSetPairs}} \code{\link{FindSynteny}}
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
                                PIDs = TRUE,
                                AcceptContigNames = TRUE,
                                Verbose = TRUE)
NewPairs <- BlockExpansion(PairsObject,
                           GapTolerance = 2L,
                           DropSingletons = FALSE,
                           Floor = 0.5,
                           Verbose = FALSE,
                           NewPairsOnly = TRUE,
                           DBPATH = DBPATH)
}
}
