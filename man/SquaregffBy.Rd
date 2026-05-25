\name{SquaregffBy}
\alias{SquaregffBy}
\title{
  Convert a GRanges object to a DataFrame
}
\description{
  A DataFrame representing a portion of the data contained in a gff file / GRanges object.
}
\usage{
  SquaregffBy(gff_object,
              collect_by = c("gene",
                             "pseudogene",
                             "ncRNA_gene",
                             "tRNA_gene"),
              verbose = FALSE)
}
\arguments{
  \item{gff_object}{
    A \code{GRanges} object produced by \code{rtracklayer}'s \code{import}.
  }
  \item{collect_by}{
    A character vector containing feature types to collect from the gff object.
  }
  \item{verbose}{
    A logical indicating whether to print out messages and a progress bar.
  }
}
\details{
  Given a DataFrame with a \code{Range} column containing an \code{IRangesList}, return a data.frame the length of the unlisted \code{IRangesList}. Each row represents the coordinates of a subfeature as described by the contents of the \code{IRangesList}. Key and SubKey columns can be used to trace child features back to their parents. 
}
\value{
  A DataFrame containing excerpts from the attributes fields of a gff file, organized by the parent feature.
}
\author{
  Nicholas Cooley \email{Nicholas.Cooley@ul.ie}
}

\seealso{
  \itemize{
    \item \code{\link{SquaregffBy}}
    \item \code{\link{FindSynteny}}
    \item \code{\link{NucleotideOverlap}}
  }
}
\examples{
library(rtracklayer)
grange_obj <- import(con = system.file("extdata",
                                       "GCF_023585725.1_ASM2358572v1_genomic.gff.gz",
                                       package = "SynExtend"),
                     format = "gff")
ImportedGFF <- SquaregffBy(gff_object = grange_obj,
                           verbose = TRUE)
}
