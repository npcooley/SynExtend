\name{ProtWeaver-ResiduePreds}
\alias{ProtWeaver-ResiduePreds}
\alias{ResidueMI.ProtWeaver}
\docType{data}
\title{
Residue Level Predictions for ProtWeaver
}
\description{
\code{ProtWeaver} incorporates four classes of prediction, each with multiple 
methods and algorithms. Residue level methods examine conservation of 
individual base pairs, commonly exhibited due to physical interactions
between proteins.

\code{predict.ProtWeaver} currently supports two Residue methods:
\itemize{
  \item \code{'ResidueMI'}
  \item \code{'NVDT'}
  \item \code{'Ancestral'}
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
base pair between two genomes. \code{ResidueMI} implements an optimization
detailed in Gao et al. (2018) to compress sequences prior to analysis 
based on correlations. This removes base pairs unlikely to contribute 
significantly to the overall MI score.

The \code{NVDT} method uses the natural vector encoding method introduced in Zhao et al. (2022). This encodes each gene sequences as a 92-dimensional vector, with the following entries:

\deqn{N(S) = (n_A,n_C,n_G,n_T,\\
              \qquad\qquad\;\,\mu_A,\mu_C,\mu_G,\mu_T,\\
              \qquad\qquad\quad\, D_2^A,D_2^C,D_2^G,D_2^T,\\
              \qquad\qquad\qquad n_{AA},n_{AC},\dots,n_{TT},\\
              \qquad\qquad\qquad\quad\;\; n_{AAA},n_{AAC},\dots,n_{TTT})}

Here \eqn{n_X} is the raw total count of nucleotide \eqn{X} (or di/trinucleotide). For single nucleotides, we also calculate \eqn{\mu_X}, the mean location of nucleotide \eqn{X}, and \eqn{D_2^X}, the second moment of the location of nucleotide \eqn{X}. The overall natural vector for a COG is calculated as the normalized mean vector from the natural vectors of all component gene sequences. Interaction scores are computed using Pearson's R between each COG's natural vector.

The \code{Ancestral} method calculates coevolution by looking at correlation of residue mutations near the leaves of each respective gene tree. 
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\references{
Gao, C. Y., et al., \emph{Correlation-compressed direct-coupling analysis}. Physical Review
E, 2018. \bold{98}(3): 032407.

Pazos, F. and A. Valencia, \emph{In silico two‐hybrid system for the selection of physically interacting protein pairs}. Proteins: Structure, Function, and Bioinformatics, 2002. \bold{47}(2): p. 219-227.

Zhao, N., et al., \emph{Protein-protein interaction and non-interaction predictions using gene sequence natural vector}. Nature Communications Biology, 2022. \bold{5}(652).
}

\seealso{

\code{\link{ProtWeaver}}

\code{\link{predict.ProtWeaver}}

\link[=ProtWeaver-PAPreds]{ProtWeaver Presence/Absence Predictors}

\link[=ProtWeaver-DMPreds]{ProtWeaver Distance Matrix Predictors}

\link[=ProtWeaver-ColocPreds]{ProtWeaver Co-localization Predictors}
}