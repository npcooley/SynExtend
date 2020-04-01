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
                  LimitIndex = FALSE,
                  OutputFormat = "Normal",
                  Verbose = FALSE)
}
\arguments{
  \item{SyntenyObject}{
An object of class ''Synteny''
}
  \item{GeneCalls}{
A list of DataFrames that describe feature calls for the sequences represented in \code{SyntenyObject}. This list must be named, and these names should match the identifiers of the associated sequences in the DECIPHER sqlite database used to build the associated synteny object. Dataframes built by ''gffToDataFrame'' can be used directly, while ''GRanges'' objects may also be used with limited functionality. Using a ''GRanges'' object will force \code{LimitIndex} to \code{FALSE}. Additionally, genomic features that extend past the terminal position of the origin FASTA sequence are truncated.
}
  \item{LimitIndex}{
Logical indicating whether to only limit which indices in a synteny object to query. \code{FALSE} by default, when \code{TRUE} only the first sequence in all selected identifiers will be used. \code{LimitIndex} can be used to skip analysis of plasmids, or solely query a single chromosome.
}
  \item{OutputFormat}{
Character string to designate how much information to return. "Sparse" returns only a filled upper triangle of exactly matched positions. "Normal" returns a matrix with associated match information in both the upper and lower triangle of the returned matrix, while "Comprehensive" will return \code{GeneCalls} used in construction in the diagonal.
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
None.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
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
