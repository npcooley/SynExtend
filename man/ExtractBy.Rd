\name{ExtractBy}
\alias{ExtractBy}
\title{
Extract and organize \code{XStringSet}s of sequences represented in a \code{PairSummaries} object.
}
\description{
Takes in a \code{PairSummaries} object and an optional vector of cluster representatives. Return an \code{XStringSet} of the sequences present in the \code{PairSummaries}, or when cluster representatives are provided, a list of \code{XStringSet}s of the sequences that make up the provided clusters.
}
\usage{
ExtractBy(x,
          y = NULL,
          DBPATH,
          Method = "all",
          DefaultTranslationTable = "11",
          Translate = TRUE,
          Storage = 1,
          Verbose = FALSE)
}
\arguments{
  \item{x}{
A \code{PairSummaries} object.
}
  \item{y}{
An optional \code{list} containing the ids of sequences in the \code{PairSummaries} object.
}
  \item{DBPATH}{
A SQLite connection object or a character string specifying the path to the database file. Constructed from DECIPHER's \code{Seqs2DB} function.
}
  \item{Method}{
How to extract sequences from the \code{PairSummaries} object. Currently only the methods ``all'' and ``clusters'' are supported.
}
  \item{Translate}{
If \code{TRUE} return \code{AAStringSet}s where possible.
}
  \item{DefaultTranslationTable}{
Currently Not Implemented! When implemented will allow for designation of a specific translation table if one is not indicated in the \code{GeneCalls} attribute of the \code{PairSummaries} object.
}
  \item{Storage}{
Numeric indicating the approximate size a user wishes to allow for holding \code{StringSet}s in memory to extract gene sequences, in ``Gigabytes''. The lower \code{Storage} is set, the more likely that \code{ExtractBy} will need to reaccess \code{StringSet}s when extracting gene sequences. The higher \code{Storage} is set, the more sequences \code{ExtractBy} will attempt to hold in memory, avoiding the need to re-access the source database many times. Set to 1 by default, indicating that \code{ExtractBy} can store a ``Gigabyte'' of sequences in memory at a time.
}
  \item{Verbose}{
Logical indicating whether to print progress bars and messages. Defaults to \code{FALSE}.
}
}
\details{
Takes in a \code{PairSummaries} object and an optional vector of cluster representatives. Return an \code{XStringSet} of the sequences present in the \code{PairSummaries}, or when cluster representatives are provided, a list of \code{XStringSet}s of the sequences that make up the provided clusters.
}
\value{
Returns either a \code{XStringSet} or a list of \code{XStringSet}s.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{FindSynteny}}, \code{\link{Synteny-class}}, \code{\link{PairSummaries}}, \code{\link{DisjointSet}}
}
\examples{
DBPATH <- system.file("extdata",
                      "Nitrosocosmicus.sqlite",
                      package = "SynExtend")
data("Nitrosocosmicus_Pairs03", package = "SynExtend")
data("Nitrosocosmicus_Sets", package = "SynExtend")

# extract the first 10 disjoint sets
Sets <- ExtractBy(x = Nitrosocosmicus_Pairs03,
                  y = Nitrosocosmicus_Sets[1:10],
                  Method = "clusters",
                  DBPATH = DBPATH,
                  Verbose = TRUE)
}
