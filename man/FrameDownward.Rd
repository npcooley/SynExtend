\name{FrameDownward}
\alias{FrameDownward}
\title{
  Adjust a DataFrame of genecalls
}
\description{
  Reframe a DataFrame of genecalls into an explicitly square representation of the features.
}
\usage{
  FrameDownward(genecalls)
}
\arguments{
  \item{genecalls}{
    A DataFrame object created by \code{\link{SquaregffBy}}.
  }
}
\details{
  Given a DataFrame with a \code{Range} column containing an \code{IRangesList}, return a data.frame the length of the unlisted \code{IRangesList}. Each row represents the coordinates of a subfeature as described by the contents of the \code{IRangesList}. Key and SubKey columns can be used to trace child features back to their parents. 
}
\value{
  A data.frame of of feature positions.
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
res <- FrameDownward(ImportedGFF)
}
