\name{SummarizePairs}
\alias{SummarizePairs}
\title{
Provide summaries of hypothetical orthologs.
}
\description{
Given the correct set of \code{SynExtend} objects and a \code{DECIPHER} database, return a data.frame of summarized genomic feature pairs. \code{SummarizePairs} will collect all the linked genomic features in the supplied \code{\link{LinkedPairs-class}} object and return descriptions of the alignments of those features.
}
\usage{
SummarizePairs(SynExtendObject,
               FeatureSeqs,
               DataBase,
               AlignmentFun = "AlignProfiles",
               RetainAnchors = FALSE,
               DefaultTranslationTable = "11",
               KmerSize = 5,
               IgnoreDefaultStringSet = FALSE,
               Verbose = FALSE,
               ShowPlot = FALSE,
               Processors = 1,
               ...)
}
\arguments{
  \item{SynExtendObject}{
An object of class \code{LinkedPairs-class}.
}
  \item{FeatureSeqs}{
An object of class \code{FeatureSeqs}.
}
  \item{DataBase}{
A character string pointing to a SQLite database, or a connection to a \code{DECIPHER} database.
}
  \item{AlignmentFun}{
A character string specifying a \code{link{DECIPHER}} alignment function. Currently only supports \code{\link{AlignProfiles}} and \code{\link{AlignPairs}}.
}
  \item{RetainAnchors}{
An argument that only affects \code{\link{AlignPairs}}; provide the kmer hits supplied by \code{\link{FindSynteny}} as alignment anchors.
}
  \item{DefaultTranslationTable}{
A character vector of length 1 identifying the translation table to use if one is not supplied in the \code{GeneCalls} attribute.
}
  \item{KmerSize}{
An integer specifying what Kmer size to collect Kmer distance between sequences at.
}
  \item{IgnoreDefaultStringSet}{
A soft memory limit for how much space to allow when building the resulting object. Translated to Gb.
}
  \item{Verbose}{
Logical indicating whether or not to display a progress bar and print the time difference upon completion.
}
  \item{ShowPlot}{
Logical indicating whether or not to provide a plot of features collected by the function.
}
  \item{Processors}{
An integer value indicating how many processors to supply to \code{\link{AlignPairs}}.
}
  \item{...}{
Additional arguments to pass to interior functions. Currently not implemented.
}
}
\details{
\code{SummarizePairs} collects features describing each linked feature pair. These include an alignment PID, an alignment Score, a Kmer distance, a concensus score for the linking hits --or whether or not linking hits are in similar places in each feature-- and a few other features.
}
\value{
An object of class \code{PairSummaries}.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{PrepareSeqs}}, \code{\link{NucleotideOverlap}}, \code{\link{FindSynteny}}, \code{\link{LinkedPairs-class}}
}
\examples{
DBPATH <- system.file("extdata",
                      "Endosymbionts_v02.sqlite",
                      package = "SynExtend")
                      
data("Endosymbionts_LinkedFeatures", package = "SynExtend")
Endosymbiont_Seqs <- PrepareSeqs(SynExtendObject = Endosymbionts_LinkedFeatures,
                                 DataBase = DBPATH,
                                 Verbose = TRUE)
SummarizedPairs <- SummarizePairs(SynExtendObject = Endosymbionts_LinkedFeatures,
                                  FeatureSeqs = Endosymbiont_Seqs,
                                  DataBase = DBPATH)
                           
}
