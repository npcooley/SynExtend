\name{CreateDecoys}
\alias{CreateDecoys}
\title{
Generating Decoy Alignments for Feature Pairs
}
\description{
A function for generating a set of decoy pairwise alignments from a database of genome sequences and an accompanying gene calls object. For each pair of genomes, one feature set is reversed by default prior to alignment, producing spurious matches that serve as a null distribution against which real \code{PairSummaries} results can be calibrated or filtered.
}
\usage{
CreateDecoys(DataBase01,
             GeneCalls,
             K_val_01,
             K_val_02,
             K_val_03,
             DefaultTranslationTable = "11",
             NT_limit = NULL,
             AA_limit = NULL,
             TRY_limit = NULL,
             DecoyMethod = c("reverse", "none"),
             Verbose = FALSE)
}
\arguments{
  \item{DataBase01}{
Either a connection object to a \code{DECIPHER}-compatible SQLite database, or a character string giving the path to such a database on disk. The database must contain a \code{Seqs} table with nucleotide sequences whose identifiers correspond to the names of \code{GeneCalls}. If amino acid sequences for a given identifier are absent will be translated on the fly and written back to the database. If a path string is supplied, the \code{RSQLite} package must be installed.
}
  \item{GeneCalls}{
A named list of gene call objects, where each name is an identifier matching an entry in the \code{Seqs} table of \code{DataBase01}. Must contain at least two entries. Each element is expected to have the fields \code{Range}, \code{Coding}, \code{Strand}, and \code{Translation_Table}, as produced by \code{SquaregffBy} or \code{FindGenes}.
}
  \item{K_val_01}{
A positive integer specifying the k-mer width used to compute nucleotide frequency profiles (\code{oligonucleotideFrequency}) for k-mer distance calculations. This value is also stored as the \code{KmerSize} attribute of the returned object.
}
  \item{K_val_02}{
A positive integer specifying the k-mer width passed to \code{IndexSeqs} when building the inverted index over amino acid feature sequences. Controls the sensitivity of the amino acid search step.
}
  \item{K_val_03}{
A positive integer specifying the k-mer width passed to \code{IndexSeqs} when building the inverted index over nucleotide feature sequences. Controls the sensitivity of the nucleotide search step.
}
  \item{DefaultTranslationTable}{
A character string of length 1 specifying the NCBI genetic code identifier to use when translating coding sequences whose translation table is not recorded in \code{GeneCalls}. Defaults to \code{"11"} (the bacterial, archaeal, and plant plastid code). Must be a valid identifier accepted by \code{Biostrings::getGeneticCode}.
}
  \item{NT_limit}{
Either \code{NULL} (the default) or a positive integer specifying the cumulative number of nucleotide-aligned decoy pairs at which iteration stops early. When \code{NULL}, no nucleotide-count limit is applied. When used concurrently with \code{AA_limit} the first limit reached triggers the early exit, leaving the second un-respected.
}
  \item{AA_limit}{
Either \code{NULL} (the default) or a positive integer specifying the cumulative number of amino-acid-aligned decoy pairs at which iteration stops early. When \code{NULL}, no amino acid count limit is applied. When used concurrently with \code{NT_limit} the first limit reached triggers the early exit, leaving the second un-respected.
}
  \item{TRY_limit}{
Either \code{NULL} (the default) or a numeric or integer of length 1 specifying the maximum number of genome-pair comparisons to attempt. When \code{NULL}, all pairs in the upper triangle of \code{GeneCalls} are attempted. Useful for generating a fixed-size decoy set from a large multi-genome collection without iterating over every possible pair.
}
  \item{DecoyMethod}{
A character string controlling how the decoy is constructed for each genome pair. Must be one of \code{"reverse"} (the default; the feature sequences of a randomly chosen genome in each pair are reversed before alignment, producing alignments expected to be spurious) or \code{"none"} (no transformation is applied, yielding alignments between unmodified sequences). Partial matching is supported.
}
  \item{Verbose}{
Logical indicating whether to print a status message, display a progress bar, and print the elapsed time upon completion. Defaults to \code{FALSE}.
}
}
\details{
For each ordered pair of genomes \eqn{(i, j)} in the upper triangle of \code{GeneCalls} (or up to \code{TRY_limit} pairs), the function:

\enumerate{
  \item Retrieves nucleotide feature sequences from \code{DataBase01} via \code{FeaturesFromDF}. Amino acid sequences are retrieved from the \code{AAs} table if present, or translated via \code{GetTranslatedFeatures} and written back to the database for reuse.
  \item Randomly assigns which genome's sequences are treated as the \dQuote{subject} and which as the \dQuote{pattern}. Depending on \code{DecoyMethod}, the subject sequences are reversed (character-level reversal of the \code{XStringSet}) before search and alignment, ensuring that any matches found are unlikely to reflect true homology.
  \item Builds inverted k-mer indices over the subject nucleotide and amino acid sequences using \code{IndexSeqs} with \code{K_val_03} and \code{K_val_02} respectively, then searches the pattern sequences against those indices using \code{SearchIndex}.
  \item For each set of search hits, computes k-mer distances (normalised Euclidean distance between \code{K_val_01}-mer frequency profiles), hit statistics (\code{TotalMatch}, \code{MaxMatch}, \code{UniqueMatches}), pairwise alignments via \code{AlignPairs}, local and approximate global PID and alignment score, a positional consensus score, and a background-corrected score using per-genome residue composition and the PFASUM50 (amino acid) or a \eqn{2/-1/1} match/mismatch/gap (nucleotide) substitution matrix.
  \item Nucleotide-level results for any pair also represented in the amino acid results are discarded to avoid redundancy.
  \item Iteration halts early if \code{AA_limit}, \code{NT_limit}, or \code{TRY_limit} is reached, whichever comes first.
}

The returned object has the same column schema as \code{PairSummaries} objects produced by \code{SummarizePairs}, making the two directly comparable for score distribution analysis and threshold calibration. All \code{Block_UID} values are assigned as sequential integers; no syntenic block structure is inferred because decoy pairs have no context from the perspective of the search scheme.

Note that this function does not implement the memory-pool management present in \code{SummarizePairs} and \code{NucleotideOverlap}, and may therefore consume more memory when applied to large genome collections.
}
\value{
A \code{data.frame} of class \code{c("data.frame", "PairSummaries")} with one row per decoy feature pair. Column definitions match those returned by \code{SummarizePairs}:

\describe{
  \item{p1}{Character. Name of the subject feature (potentially reversed, depending on \code{DecoyMethod}).}
  \item{p2}{Character. Name of the pattern feature.}
  \item{Consensus}{Numeric. Mean positional consensus score across k-mer hits for the pair, re-scaled so that values near 1 indicate hits at similar relative positions in both sequences and values near 0 indicate maximal positional disagreement.}
  \item{p1featurelength}{Integer. Length of the subject feature in nucleotides.}
  \item{p2featurelength}{Integer. Length of the pattern feature in nucleotides.}
  \item{blocksize}{Integer. Always 1 for decoy pairs; no syntenic block context is available.}
  \item{KDist}{Numeric. Normalised Euclidean distance between the \code{K_val_01}-mer frequency profiles of the two features.}
  \item{TotalMatch}{Integer. Total number of nucleotide positions covered by k-mer hits for the pair.}
  \item{MaxMatch}{Integer. Length of the largest single k-mer hit for the pair, in nucleotides.}
  \item{UniqueMatches}{Integer. Number of distinct k-mer hits for the pair.}
  \item{Local_PID}{Numeric. Fraction of matched positions within the aligned region (local percent identity).}
  \item{Local_Score}{Numeric. Alignment score normalised by alignment length.}
  \item{Approx_Global_PID}{Numeric. Fraction of matched positions normalised by the length of the longer feature (approximate global percent identity).}
  \item{Approx_Global_Score}{Numeric. Alignment score normalised by the length of the longer feature.}
  \item{Alignment}{Character. Alphabet used for alignment: \code{"AA"} for amino acid or \code{"NT"} for nucleotide.}
  \item{Block_UID}{Integer. A sequential unique identifier assigned to each row; carries no syntenic meaning for decoy pairs.}
  \item{Delta_Background}{Numeric. The difference between the approximate global alignment score and the expected score under a background model derived from genome-wide residue composition.}
}

The returned object retains the following attributes:
\describe{
  \item{GeneCalls}{The \code{GeneCalls} list supplied to the function.}
  \item{KmerSize}{The value of \code{K_val_01} used during the run.}
}
}
\author{
  Nicholas Cooley \email{Nicholas.Cooley@ul.ie}
}
\seealso{
  \itemize{
    \item \code{\link{SummarizePairs}}
    \item \code{\link{NucleotideOverlap}}
    \item \code{\link{EvaluatePairs}}
  }
}
\examples{
library(DBI)
data("genecalls")
tmp01 <- system.file("extdata",
                      "example_db.sqlite",
                      package = "SynExtend")
tmp02 <- tempfile()
file.copy(from = tmp01,
          to = tmp02)

drv <- dbDriver("SQLite")
conn01 <- dbConnect(drv = drv,
                    tmp02)

x <- CreateDecoys(DataBase01 = conn01,
                  GeneCalls = genecalls,
                  K_val_01 = 5,
                  K_val_02 = 5,
                  K_val_03 = 10,
                  DefaultTranslationTable = "11",
                  TRY_limit = 1,
                  DecoyMethod = c("reverse", "none"),
                  Verbose = FALSE)
}
