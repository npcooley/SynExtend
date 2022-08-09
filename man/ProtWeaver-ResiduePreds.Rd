\name{ProtWeaver-ResiduePreds}
\alias{ProtWeaver-ResiduePreds}
\alias{ResidueMI.ProtWeaver}
\docType{data}
\title{
Residue Based Predictions for ProtWeaver
}
\description{
\code{ProtWeaver} incorporates four classes of prediction, each with multiple 
methods and algorithms. Residue level methods examine conservation of 
individual base pairs, commonly exhibited due to physical interactions
between proteins.

\code{predict.ProtWeaver} currently supports one Residue method:
\itemize{
  \item \code{'ResidueMI'}
}
}

\details{
All residue methods require a \code{ProtWeaver} object initialized
with \code{dendrogram} objects and ancestral states. See \code{\link{ProtWeaver}} for more information on input data types.

The \code{ResidueMI} method looks at mutual information between 
ancestral gene states at the nucleotide level. This approach is conceptually 
similar to the \emph{in silico} two-hybrid model introduced in Pazos and 
Valencia (2002). Gene ancestral states are paired based on how they partition
the tree; the resulting states are concatenated and analyzed. 

Residue level methods scale poorly because they need to analyze every
base pair between two genomes. \code{ProtWeaver} implements an optimization
detailed in Gao et al. (2018) to compress sequences prior to analysis 
based on correlations. This removes base pairs unlikely to contribute 
significantly to the overall MI score.
}

\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\references{
Gao, C. Y., et al., \emph{Correlation-compressed direct-coupling analysis}. Physical Review
E, 2018. \bold{98}(3): 032407.

Pazos, F. and A. Valencia, \emph{In silico two‐hybrid system for the selection of physically interacting protein pairs}. Proteins: Structure, Function, and Bioinformatics, 2002. \bold{47}(2): p. 219-227.
}

\seealso{

\code{\link{ProtWeaver}}

\code{\link{predict.ProtWeaver}}

\link[=ProtWeaver-PAPreds]{ProtWeaver Presence/Absence Predictors}

\link[=ProtWeaver-DMPreds]{ProtWeaver Distance Matrix Predictors}

\link[=ProtWeaver-ColocPreds]{ProtWeaver Co-localization Predictors}
}