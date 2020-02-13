\name{gffToDataFrame}
\alias{gffToDataFrame}
\title{
Generate a data.frame of gene calls from a gff3 file
}
\description{
Generate a data.frame of gene calls from a gff3 file
}
\usage{
gffToDataFrame(GFF,
               AdditionalAttrs = NULL,
               AdditionalTypes = NULL,
               RawTableOnly = FALSE,
               Verbose = FALSE)
}

\arguments{
  \item{GFF}{
A url or filepath specifying a gff3 file to import
}
  \item{AdditionalAttrs}{
A vector of character strings to expand the attributes to pull.
}
  \item{AdditionalTypes}{
A vector of character strings to query from the the ''Types'' column.
}
  \item{RawTableOnly}{
Logical specifying whether to return the raw imported GFF without complex parsing.
}
  \item{Verbose}{
Logical specifying whether to print a progress bar and time difference.
}
}
\details{
Import a gff file into a rectangular parsable object.
}
\value{
A data.frame
}
\references{
None yet!
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}
\note{
None.
}

\seealso{

}
\examples{
ImportedGFF <- gffToDataFrame(GFF = system.file("extdata",
                                                "GCA_000269985.1_ASM26998v1_genomic.gff.gz",
                                                package = "SynExtend"),
                              Verbose = TRUE)
}

\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
