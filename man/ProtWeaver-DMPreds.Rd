\name{ProtWeaver-DMPreds}
\alias{ProtWeaver-DMPreds}
\alias{MirrorTree.ProtWeaver}
\alias{ContextTree.ProtWeaver}
\docType{data}
\title{
Distance Matrix Predictions for ProtWeaver
}
\description{
\code{ProtWeaver} incorporates four classes of prediction, each with multiple 
methods and algorithms. Distance Matrix (DM) methods examine conservation 
of overall evolutionary rates within orthology groups using distance matrices
constructed from each gene tree.

\code{predict.ProtWeaver} currently supports three DM methods:
\itemize{
  \item \code{'MirrorTree'}
  \item \code{'ContextTree'}
  \item \code{'TreeDistance'}
}
}

\details{
All distance matrix methods require a \code{ProtWeaver} object initialized
with \code{dendrogram} objects. See \code{\link{ProtWeaver}} for more information on input data types.

The \code{MirrorTree} method was introduced by Pazos et al. (2001). This method
builds distance matrices using a nucleotide substitution model, and then
calculates coevolution between gene families using the Pearson correlation 
coefficient of the upper triangle of the two corresponding matrices.

Subsequent work by Pazos et al. (2005) and Sato et al. (2005, 2006) found
multiple ways to improve predictions from the initial \code{MirrorTree} method.
These methods incorporate additional phylogenetic context, and are thus
called \code{ContextTree} methods. These improvements include correcting
for overall evolutionary rate using a species tree or projection vectors, using partial correlation coefficient instead of Pearson correlation, and/or 
normalizing distances prior to calculation. The built-in \code{ContextTree}
method uses corrections based on what is possible with provided input data.

The \code{TreeDistance} method uses phylogenetic tree distance to quantify differences between gene trees. This method implements a number of metrics and groups them together to improve overall runtime. Individual methods can be specified using the \code{TreeMethods} argument, which expects a character vector containing one or more of the following:
\itemize{
\item \code{"GRF"}: \link[=GRFDist]{Information-Theoretic Robinson-Foulds Distance}
\item \code{"RF"}: \link[=RFDist]{Robinson-Foulds Distance}
\item \code{"JRF"}: \link[=JRFDist]{Jaccard-Robinson-Foulds Distance}
\item \code{"Nye"}: \link[=JRFDist]{Nye Similarity}
\item \code{"KF"}: \link[=KFDist]{Kuhner-Felsenstein Distance}
\item \code{"all"}: All of the above methods
}

See the links above for more information and references. All of these metrics are accessible using the \code{\link{PhyloDistance}} method. Method \code{"JRF"} defaults to a \code{k}
value of 2, but this can be specified further if necessary using the \code{JRFk} input parameter. Higher values of \code{k} approach the value of Robinson-Foulds 
distance, but these have a negligible impact on performance so use of the default parameter is encouraged for simplicity.
}

\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\references{
Pazos, F. and A. Valencia, \emph{Similarity of phylogenetic trees as indicator of protein–protein interaction.} Protein Engineering, Design and Selection, 2001. \bold{14}(9): p. 609-614.

Pazos, F., et al., \emph{Assessing protein co-evolution in the context of the tree of life assists in the prediction of the interactome.} J Mol Biol, 2005. \bold{352}(4): p. 1002-15.

Sato, T., et al., \emph{The inference of protein-protein interactions by co-evolutionary analysis is improved by excluding the information about the phylogenetic relationships.} Bioinformatics, 2005. \bold{21}(17): p. 3482-9.

Sato, T., et al., \emph{Partial correlation coefficient between distance matrices as a new indicator of protein-protein interactions.} Bioinformatics, 2006. \bold{22}(20): p. 2488-92.
}

\seealso{

\code{\link{ProtWeaver}}

\code{\link{predict.ProtWeaver}}

\link[=ProtWeaver-PAPreds]{ProtWeaver Presence/Absence Predictors}

\link[=ProtWeaver-ColocPreds]{ProtWeaver Co-localization Predictors}

\link[=ProtWeaver-ResiduePreds]{ProtWeaver Residue Level Predictors}

\code{\link{PhyloDistance}}
}