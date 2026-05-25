\name{GetTranslatedFeatures}
\alias{GetTranslatedFeatures}
\title{
Translate Nucleotide Features to Amino Acid Sequences
}
\description{
A lightweight function for translating a set of nucleotide feature sequences into amino acid sequences. Only features that are flagged as coding and whose total coding length is in phase will be translated. Multiple genetic codes may be applied in a single call when different features carry different translation table annotations. The original feature order is preserved in the returned \code{AAStringSet}.
}
\usage{
GetTranslatedFeatures(Nucs,
                      GeneCalls,
                      DefaultTranslationTable = "11")
}
\arguments{
  \item{Nucs}{
A \code{DNAStringSet} of nucleotide feature sequences to be translated, typically produced by \code{FeaturesFromDF}. Names are expected to follow the internal \code{SynExtend} naming convention (underscore-delimited fields, with the feature index in the third position) so that output ordering can be restored after grouping by genetic code.
}
  \item{GeneCalls}{
A \code{DataFrame} of gene calls corresponding to the features in \code{Nucs}. Must contain at least the following columns:
\describe{
  \item{\code{Coding}}{Logical vector indicating whether each feature is a coding sequence.}
  \item{\code{Translation_Table}}{Character vector of NCBI genetic code identifiers (e.g. \code{"11"}), one per feature. \code{NA} values are replaced by \code{DefaultTranslationTable} for features that are both coding and in-frame.}
  \item{\code{Range}}{A list of \code{IRanges} or \code{GRanges} objects giving the genomic coordinates of each feature. The sum of widths for each element is used to determine whether the feature length is a multiple of three.}
}
}
  \item{DefaultTranslationTable}{
A character string of length 1 specifying the NCBI genetic code identifier to use for coding, in-frame features whose \code{Translation_Table} entry is \code{NA}. Defaults to \code{"11"} (the bacterial, archaeal, and plant plastid code). Must be a valid identifier accepted by \code{Biostrings::getGeneticCode}.
}
}
\details{
Translation is restricted to features satisfying all three conditions: the feature is flagged as coding (\code{GeneCalls$Coding == TRUE}), its total nucleotide length is a multiple of three (sum of \code{width(GeneCalls$Range)} modulo 3 equals zero), and a translation table is available after substituting \code{DefaultTranslationTable} for any \code{NA} entries.

Features meeting these criteria are grouped by their genetic code identifier and translated together using \code{Biostrings::translate} with \code{if.fuzzy.codon = "solve"}. When more than one genetic code is present, the per-code results are concatenated and reordered to match the original position order of \code{Nucs}, relying on the feature index encoded in the third underscore-delimited field of each sequence name.

Features that are non-coding or not in-frame are silently excluded; the returned \code{AAStringSet} therefore contains fewer sequences than \code{Nucs} whenever such features are present.
}
\value{
An \code{AAStringSet} containing the translated amino acid sequences for all coding, in-frame features in \code{Nucs}, in the same relative order as their source sequences appear in \code{Nucs}. Names are inherited from the corresponding entries in \code{Nucs}.
}
\author{
  Nicholas Cooley \email{Nicholas.Cooley@ul.ie}
}
\seealso{
  \itemize{
    \item \code{\link{FeaturesFromDF}}
    \item \code{\link{SummarizePairs}}
    \item \code{\link{CreateDecoys}}
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
x <- SearchDB(dbFile = conn01,
              identifier = "1",
              nameBy = "description")

y <- FeaturesFromDF(Genome = x,
                    GeneCalls = genecalls[[1]],
                    Index = "1")
                    
z <- GetTranslatedFeatures(Nucs = y,
                           GeneCalls = genecalls[[1]],
                           DefaultTranslationTable = "11")
}
