\name{PrepareSeqs}
\alias{PrepareSeqs}
\title{
Add feature sequences to Decipher databases.
}
\description{
Given a \code{SynExtend} object with a \code{GeneCalls} attribute, and a \code{DECIPHER} database, add sequence tables named 'AAs' and 'NTs' to the database. The new tables contain all translatable sequences indicated by the genecalls, and all nucleotide feature sequences.
}
\usage{
PrepareSeqs(SynExtendObject,
            DataBase01,
            DefaultTranslationTable = "11",
            Identifiers = NULL,
            Verbose = FALSE)
}
\arguments{
  \item{SynExtendObject}{
An object of class \code{PairSummaries} or of \code{LinkedPairs}. Object must have a \code{GeneCalls} attribute.
}
  \item{DataBase01}{
A character string pointing to a SQLite database, or a connection to a \code{DECIPHER} database.
}
  \item{DefaultTranslationTable}{
A character vector of length 1 identifying the translation table to use if one is not supplied in the \code{GeneCalls} attribute.
}
  \item{Identifiers}{
By default \code{NULL}, but can be used to supply a vector of character identifiers for returning a subset of prepared sequences.
}
  \item{Verbose}{
Logical indicating whether or not to display a progress bar and print the time difference upon completion.
}
}
\details{
\code{PrepareSeqs} adds two tables to a DECIPHER database. One named 'AAs' that contains all translatable features, i.e. features with a coding length divisible by 3 and designated as coding. And another named 'NTs' which contains all features.
}
\value{
An integer count of the number of feature sets added to the DECIPHER database.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{SummarizePairs}}, \code{\link{NucleotideOverlap}}, \code{\link{FindSynteny}}
}
\examples{
DBPATH <- system.file("extdata",
                      "Endosymbionts_v02.sqlite",
                      package = "SynExtend")
                      
data("Endosymbionts_LinkedFeatures", package = "SynExtend")
# this will add seqs to the DB
# PrepareSeqs(SynExtendObject = Endosymbionts_LinkedFeatures,
#             DataBase = DBPATH,
#             Verbose = TRUE)
}
