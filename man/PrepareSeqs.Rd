\name{PrepareSeqs}
\alias{PrepareSeqs}
\title{
Return gene sequences.
}
\description{
Given a \code{SynExtend} object with a \code{GeneCalls} attribute, and a \code{DECIPHER} database, return all gene sequences and their translations.
}
\usage{
PrepareSeqs(SynExtendObject,
            DataBase,
            DefaultTranslationTable = "11",
            Identifiers = NULL,
            Storage = 1,
            Verbose = FALSE)
}
\arguments{
  \item{SynExtendObject}{
An object of class \code{PairSummaries} or of \code{LinkedPairs}. Object must have a \code{GeneCalls} attribute.
}
  \item{DataBase}{
A character string pointing to a SQLite database, or a connection to a \code{DECIPHER} database.
}
  \item{DefaultTranslationTable}{
A character vector of length 1 identifying the translation table to use if one is not supplied in the \code{GeneCalls} attribute.
}
  \item{Identifiers}{
By default \code{NULL}, but can be used to supply a vector of character identifiers for returning a subset of prepared sequences.
}
  \item{Storage}{
A soft memory limit for how much space to allow when building the resulting object. Translated to Gb.
}
  \item{Verbose}{
Logical indicating whether or not to display a progress bar and print the time difference upon completion.
}
}
\details{
\code{PrepareSeqs} returns the sequences of genes and their translations where appropriate.
}
\value{
An object of class \code{FeatureSeqs}.
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
CurrentSeqs <- PrepareSeqs(SynExtendObject = Endosymbionts_LinkedFeatures,
                           DataBase = DBPATH,
                           Verbose = TRUE)
}
