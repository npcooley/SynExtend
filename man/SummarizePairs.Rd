\name{SummarizePairs}
\alias{SummarizePairs}
\title{
Summarizing Linked Feature Pairs from a LinkedPairs Object
}
\description{
A function for generating a comprehensive summary of linked genomic feature pairs from an object of class \code{LinkedPairs}. For each pair, the function computes alignment statistics (percent identity, alignment score), k-mer distance, syntenic block context, and an approximate background-corrected score, returning results as an object of class \code{PairSummaries}.
}
\usage{
SummarizePairs(SynExtendObject,
               DataBase01,
               DefaultTranslationTable = "11",
               KmerSize = 5,
               Verbose = FALSE,
               ShowPlot = FALSE,
               Processors = 1,
               Storage = 2,
               Anchors = c("enforce", "infer", "ignore"),
               ...)
}
\arguments{
  \item{SynExtendObject}{
An object of class \code{LinkedPairs}, typically produced by \code{NucleotideOverlap}. Each populated cell in the upper triangle of this matrix object describes the syntenic hits that link a pair of genomes.
}
  \item{DataBase01}{
Either a connection object to a \code{DECIPHER}-compatible SQLite database, or a character string giving the path to such a database on disk. The database must contain nucleotide sequences (indexed by identifier) that correspond to the genomes referenced in \code{SynExtendObject}. If a path string is provided, the \code{RSQLite} package must be installed. Amino acid sequences will be read from the database if an \code{AAs} table is present; otherwise they will be computed on the fly and written back to the database.
}
  \item{DefaultTranslationTable}{
A character string of length 1 specifying the NCBI genetic code identifier to use when translating coding sequences whose translation table is not recorded in the \code{GeneCalls} attribute of \code{SynExtendObject}. Defaults to \code{"11"} (the bacterial, archaeal, and plant plastid code). Must be a valid identifier accepted by \code{Biostrings::getGeneticCode}.
}
  \item{KmerSize}{
A positive integer specifying the k-mer width used to compute nucleotide frequency profiles for each feature. These profiles are used to derive a k-mer distance (\code{KDist}) between each pair of features. Must be less than 10. Defaults to \code{5}.
}
  \item{Verbose}{
Logical indicating whether to display a progress bar and print the elapsed time upon completion. Defaults to \code{FALSE}.
}
  \item{ShowPlot}{
Logical. Reserved for future use. Currently has no effect. Defaults to \code{FALSE}.
}
  \item{Processors}{
A positive whole-number integer, or \code{NULL}, specifying the number of processors to use for parallel operations passed to \code{AlignPairs}. When \code{NULL}, the number of available cores is detected automatically via \code{DECIPHER}. Defaults to \code{1}.
}
  \item{Storage}{
A positive numeric value specifying the maximum amount of memory (in gigabytes) that the internal data pool may occupy before cached sequence data for earlier genomes is evicted. Increase this value when analysing many large genomes to avoid redundant database queries; decrease it on memory-constrained systems. Defaults to \code{2}.
}
  \item{Anchors}{
A character string controlling how alignment anchors are constructed when calling \code{AlignPairs}. Must be one of \code{"enforce"} (force terminal anchors at the ends of each sequence, the default), \code{"infer"} (infer anchors from syntenic hit positions; not yet implemented), or \code{"ignore"} (pass no anchors). Partial matching is supported.
}
  \item{\dots}{
Additional named arguments reserved for future use. Currently not passed to any downstream function.
}
}
\details{
For each occupied cell in the upper triangle of \code{SynExtendObject}, the function:

\enumerate{
  \item Retrieves or constructs nucleotide and amino acid feature sequences from \code{DataBase01} using the \code{GeneCalls} stored as an attribute of \code{SynExtendObject}.
  \item Computes k-mer frequency profiles (\code{KmerSize}-mers) for all nucleotide features and derives a normalised Euclidean k-mer distance for each pair.
  \item Aligns each pair with \code{AlignPairs}, optionally anchored according to the \code{Anchors} argument.
  \item Computes local and approximate global percent identity (PID) and alignment score from the \code{AlignPairs} output.
  \item Computes a per-pair approximate background alignment score (\code{Delta_Background}) by subtracting an expected score derived from the residue composition of the two genomes from the observed approximate global alignment score.
  \item Assigns syntenic block identifiers (\code{Block_UID}) via \code{BlockByRank}, with a globally unique offset applied across all pairwise comparisons so that block IDs do not collide between genome pairs.
  \item Computes a hit-weighted positional consensus score (\code{Consensus}) from the syntenic hit positions stored in the lower triangle of \code{SynExtendObject}.
}

If Amino acid sequences for a given identifier are not already present they are written back to the database under the table name \code{AAs} and reused on subsequent iterations, reducing redundant computation for large multi-genome comparisons.

Memory consumption is managed by evicting the least-recently-needed genome data from the internal pool when its size exceeds \code{Storage} gigabytes. In practice this is unlikely to matter for most users unless you are pushing the bounds of storage and memory on consumer hardware.
}
\value{
A \code{data.frame} of class \code{c("data.frame", "PairSummaries")} with one row per linked feature pair. Columns are:

\describe{
  \item{p1}{Character. Name of the query feature.}
  \item{p2}{Character. Name of the subject feature.}
  \item{Consensus}{Numeric. A positional consensus score (between 0 and 1) summarising the agreement of syntenic hit anchors across the pair; 1 indicates complete disagreement and 0 indicates perfect agreement.}
  \item{p1featurelength}{Integer. Nucleotide length of the query feature.}
  \item{p2featurelength}{Integer. Nucleotide length of the subject feature.}
  \item{blocksize}{Integer. Number of feature pairs belonging to the same syntenic block as this pair.}
  \item{KDist}{Numeric. Normalised Euclidean distance between the k-mer frequency profiles of the two features.}
  \item{TotalMatch}{Integer. Total number of exact nucleotide overlap bases recorded for this pair in the \code{LinkedPairs} object.}
  \item{MaxMatch}{Integer. Maximum single k-mer hit size recorded for this pair.}
  \item{UniqueMatches}{Integer. Number of unique k-mer matches recorded for this pair.}
  \item{Local_PID}{Numeric. Fraction of matched positions within the aligned region (local percent identity).}
  \item{Local_Score}{Numeric. Alignment score normalised by the alignment length.}
  \item{Approx_Global_PID}{Numeric. Fraction of matched positions normalised by the length of the longer feature (approximate global percent identity).}
  \item{Approx_Global_Score}{Numeric. Alignment score normalised by the length of the longer feature.}
  \item{Alignment}{Character. Alphabet used for alignment: \code{"AA"} for amino acid or \code{"NT"} for nucleotide.}
  \item{Block_UID}{Integer. A globally unique block identifier grouping syntenic feature pairs that co-occur in the same syntenic block. Pairs not assigned to a multi-pair block receive a unique singleton ID.}
  \item{Delta_Background}{Numeric. The difference between the approximate global alignment score and the expected score under a background model derived from genome-wide residue composition; higher values indicate greater similarity above background.}
}

The returned object retains the following attributes:
\describe{
  \item{GeneCalls}{The named list of gene calls inherited from \code{SynExtendObject}.}
  \item{KmerSize}{The value of \code{KmerSize} used during the run.}
  \item{DefaultTranslationTable}{The value of \code{DefaultTranslationTable} used during the run.}
}
}
\author{
  Nicholas Cooley \email{Nicholas.Cooley@ul.ie}
}
\seealso{
  \itemize{
    \item \code{\link{NucleotideOverlap}}
    \item \code{\link{FindSynteny}}
  }
}
\examples{
library(DBI)
data("linked_features")
tmp01 <- system.file("extdata",
                      "example_db.sqlite",
                      package = "SynExtend")
tmp02 <- tempfile()
file.copy(from = tmp01,
          to = tmp02)

drv <- dbDriver("SQLite")
conn01 <- dbConnect(drv = drv,
                    tmp02)

x <- SummarizePairs(SynExtendObject = linked_features,
                    DataBase01 = conn01)
}
