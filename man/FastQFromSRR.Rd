\name{FastQFromSRR}
\alias{FastQFromSRR}
\title{
Get Sequencing Data from the SRA
}
\description{
Get sequencing data from the SRA. 
}
\usage{
FastQFromSRR(SRR,
             ARGS = list("--gzip" = NULL,
                         "--skip-technical" = NULL,
                         "--readids" = NULL,
                         "--read-filter" = "pass",
                         "--dumpbase" = NULL,
                         "--split-3" = NULL,
                         "--clip" = NULL),
             KEEPFILES = FALSE)
}
\arguments{
  \item{SRR}{
A character vector of length 1 representing an SRA Run Accession, such as one that would be passed to the \code{prefetch}, \code{fastq-dump}, or \code{fasterq-dump} functions in the SRAToolkit.
}
  \item{ARGS}{
A list representing key and value sets used to construct the call to \code{fastq-dump}, multi-argument values are passed to \code{paste} directly and should be structured accordingly.
}
  \item{KEEPFILES}{
Logical indicating whether or not keep the downloaded fastq files outside of the R session. If \code{TRUE}, downloaded files will be moved to R's working directory with the default names assigned by \code{fastq-dump}. If \code{FALSE} - the default, they are removed and only the list of \code{QualityScaledDNAStringSets} returned by the function are retained.
}
}
\details{
\code{FastQFromSRR} is a barebones wrapper for \code{fastq-dump}, it is set up for convenience purposes only and does not add any additional functionality. Requires a functioning installation of the SRAtoolkit.
}
\value{
A list of \code{QualityScaledDNAStringSets}. The composition of this list will be determined by \code{fastq-dump}'s splitting arguments.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}
\examples{
\donttest{
x <- "ERR10466327"
y <- FastQFromSRR(SRR = x)
}
}
