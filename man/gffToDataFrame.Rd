\name{gffToDataFrame}
\alias{gffToDataFrame}
\title{
Generate a DataFrame of gene calls from a gff3 file
}
\description{
Generate a DataFrame of gene calls from a gff3 file
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
A vector of character strings to designate the attributes to pull. Default Attributes include: ''ID'', ''Parent'', ''Name'', ''gbkey'', ''gene'', ''product'', ''protein_id'', ''gene_biotype'', and ''Note''.
}
  \item{AdditionalTypes}{
A vector of character strings to query from the the ''Types'' column. Default types are limited to ''Gene'' and ''Pseudogene'', but any possible entry for ''Type'' in a gff3 format can be added, such as ''rRNA'', or ''CRISPR_REPEAT''.
}
  \item{RawTableOnly}{
Logical specifying whether to return the raw imported GFF without complex parsing. Remains as a holdover from function construction and debugging. For simple gff3 import see \code{rtracklayer::import}.
}
  \item{Verbose}{
Logical specifying whether to print a progress bar and time difference.
}
}
\details{
Import a gff file into a rectangular parsable object.
}
\value{
A DataFrame with relevant information extracted from a GFF.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\examples{
ImportedGFF <- gffToDataFrame(GFF = system.file("extdata",
                                                "GCA_006740685.1_ASM674068v1_genomic.gff.gz",
                                                package = "SynExtend"),
                              Verbose = TRUE)
}

\keyword{GeneCalls}
