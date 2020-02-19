\name{NucleotideOverlap}
\alias{NucleotideOverlap}
\title{
Tabulating Pairs of Genomic Sequences
}
\description{
A function for concisely tabulating where genomic features are connected by syntenic hits.
}
\usage{
NucleotideOverlap(SyntenyObject,
                  GeneCalls,
                  LimitIndex = TRUE,
                  OutputFormat = "Normal",
                  Verbose = FALSE)
}
\arguments{
  \item{SyntenyObject}{
An object of class ``Synteny''
}
  \item{GeneCalls}{
A list of DataFrames that describe feature calls for the sequences represented in \code{SyntenyObject}
}
  \item{LimitIndex}{
Logical indicating whether to only limit which indices in a synteny object to query. \code{TRUE} by default, which means that only the first indices for each cell in the synteny object are queried.
}
  \item{OutputFormat}{
Character string to designate whether to output a lot, or a little information.
}
  \item{Verbose}{
Logical indicating whether or not to display a progress bar and print the time difference upon completion.
}
}
\details{
Builds a matrix of lists that contain information about linked pairs.
}
\value{
An object of class ''LinkedPairs''.
}
\references{
None yet.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}
\note{
A note.
}


\seealso{
\code{\link{FindSynteny}}, \code{\link{Synteny-class}}
}
\examples{
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
}
