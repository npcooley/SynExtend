\name{ProtWeaver-ColocPreds}
\alias{ProtWeaver-ColocPreds}
\alias{Coloc.ProtWeaver}
\docType{data}
\title{
Co-localization Predictions for ProtWeaver
}
\description{
\code{ProtWeaver} incorporates four classes of prediction, each with multiple 
methods and algorithms. Co-localization (Coloc) methods examine
conservation of relative location and transcriptional direction of genetic
regions within the genome.

\code{predict.ProtWeaver} currently supports one Coloc method:
\itemize{
  \item \code{'Coloc'}
}
}

\details{
All distance matrix methods require a \code{ProtWeaver} object initialized
with gene locations using the a three number code. See \code{\link{ProtWeaver}} for more information on input data types.

The built-in \code{Coloc} examines relative location of genes within genomes
as evidence of interaction. For a given pair of genes, the score is given by
\eqn{\sum_{G}e^{-|dI_G|}}, where \eqn{G} the set of genomes and \eqn{dI_G} the 
difference in index between the two genes in genome \eqn{G}. Using gene index
instead of number of base pairs avoids bias introduced by gene and genome length.
}

\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\seealso{

\code{\link{ProtWeaver}}

\code{\link{predict.ProtWeaver}}

\link[=ProtWeaver-PAPreds]{ProtWeaver Presence/Absence Predictors}

\link[=ProtWeaver-DMPreds]{ProtWeaver Distance Matrix Predictors}

\link[=ProtWeaver-ResiduePreds]{ProtWeaver Residue Level Predictors}
}