\name{FindColocalSets}
\alias{FindColocalSets}
\title{
  Find co-localized groups of single-copy co-occurring sets
}
\description{
  Identify groups of co-occuring sets whose members sit on the same contig and within a maximum nucleotide distance of each other in \emph{every} genome that the sets cover.
}
\usage{
FindColocalSets(genecalls,
                detected_sets,
                allowed_miss = 0,
                allow_subset = FALSE,
                complete_set = FALSE,
                max_gap = 5000)
}
\arguments{
  \item{genecalls}{
    A named \code{list} of \code{DataFrame}s, one per genome, as produced by \code{SquaregffBy} or \code{gffToDataFrame}.
    Names must be plain integer identifiers (\code{"1"}, \code{"2"}, \ldots) matching the identifiers used when sequences were imported with \code{Seqs2DB}.
    Every element must have at least one row and the columns \code{Contig}, \code{Start}, and \code{Stop}.
}
  \item{detected_sets}{
    A \code{list} of \code{character} vectors. Each vector is one set of features, given as ids of the form \code{"genome_index_row"}.
    \code{genome} is the persistent identifier, \code{index} is the sequence index within that genome, and \code{row} is the row of the feature in \code{genecalls[[genome]]}.
    Sets must be single copy: no genome may appear more than once within a set.
}
  \item{allowed_miss}{
    Reserved for future use. Warns if not \code{0}.
}
  \item{allow_subset}{
    Reserved for future use. Warns if not \code{FALSE}.
}
  \item{complete_set}{
    Logical of length 1. If \code{TRUE}, only sets that contain exactly one feature from every genome in \code{genecalls} are evaluated.
    If \code{FALSE}, sets are grouped by the exact collection of genomes they cover, and each group of two or more sets is evaluated separately.
    \code{FALSE} by default.
}
  \item{max_gap}{
    Numeric of length 1, non-negative. The largest number of nucleotides allowed between the end of one feature and the start of the next feature in the same co-localized run.
    Overlapping features have a negative gap and are always linked.
}
}
\details{
  \code{FindColocalSets} detects where co-occuring feature sets are also co-localized.
  Testing is performed against the order of the features *on the genome* making providing greater latitute for inclusion of rearrangements.
  Future expansion of functionality to identify partially co-local sets and co-local subsets of co-occuring supersets is planned, but not implemented.
}
\value{
  A \code{list} of \code{integer} vectors. Each vector is one co-localized group and holds the positions in \code{detected_sets} of the sets in that group.
  If no co-localized groups are found, or \code{detected_sets} is empty, an empty \code{list} is returned.
}
\author{
  Nicholas Cooley \email{nicholas.cooley@ul.ie}
}
\seealso{
  \code{\link{SquaregffBy}},
  \code{\link{SummarizePairs}},
  \code{\link{EvaluatePairs}},
  \code{\link{ExoLabel}}
}
\examples{

library(DBI)
data("init_pairs")
tmp01 <- system.file("extdata",
                      "example_db.sqlite",
                      package = "SynExtend")
tmp02 <- tempfile()
file.copy(from = tmp01,
          to = tmp02)

drv <- dbDriver("SQLite")
conn01 <- dbConnect(drv = drv,
                    tmp02)
                    
x <- EvaluatePairs(InputPairs = init_pairs,
                   DataBase01 = conn01)
                   
tmp_pairs <- tempfile()
write.table(x = x[, c(1, 2, 13)],
            sep = "\t",
            quote = FALSE,
            append = FALSE,
            col.names = FALSE,
            row.names = FALSE,
            file = tmp_pairs)
y <- ExoLabel(edgelistfiles = tmp_pairs,
              return_table = TRUE)
unlink(tmp_pairs)
clusts <- y$results

clusts <- tapply(X = clusts$Vertex,
                 INDEX = clusts$Cluster,
                 FUN = c)
# index names are not meaningful
names(clusts) <- NULL

# clusters *must* be single copy
single_copy <- vapply(X = clusts,
                      FUN = function(x) {
                        y <- strsplit(x = x,
                                      split = "_",
                                      fixed = TRUE)
                        y <- do.call(rbind,
                                     y)
                        nrow(y) == length(unique(y[, 1, drop = TRUE]))
                      },
                      FUN.VALUE = vector(mode = "logical",
                                         length = 1))
clusts <- clusts[single_copy]

z <- FindColocalSets(genecalls = attr(x = init_pairs,
                                      which = "GeneCalls"),
                     detected_sets = clusts,
                     max_gap = 5000)

}
