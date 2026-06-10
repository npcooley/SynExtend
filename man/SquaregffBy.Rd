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
              additional_columns = NULL,
              verbose = FALSE)
}
\arguments{
  \item{gff_object}{
    A \code{GRanges} object produced by \code{rtracklayer}'s \code{import}.
  }
  \item{collect_by}{
    A character vector containing feature types to collect from the gff object.
  }
  \item{additional_columns}{
    An optional character vector of attribute fields to pull from the \code{GRanges} object and capture as a \code{CharacterList}. \code{NULL} by default.
  }
  \item{verbose}{
    A logical indicating whether to print out messages and a progress bar.
  }
}
\details{
  Given a \code{GRanges} object return a DataFrame where child features have been condensed so that the sequential scroll through the rows is working through only parent features in the original objects. In the case of coding features with multiple splice variants, only the first is represented.
  
  Only feature types enumerated in the \code{collect_by} argument will be allowed to be parent features. Attributes designated by the \code{additional_columns} argument will be be captured as \code{CharacterList}s.
}
\value{
  A DataFrame containing excerpts from the attributes fields of a gff file, organized by the parent feature. Feature locations are included as an \code{IRangesList}.
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
