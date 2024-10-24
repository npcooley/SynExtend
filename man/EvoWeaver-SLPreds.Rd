\name{EvoWeaver-SLPreds}
\alias{EvoWeaver-SLPreds}
\alias{SequenceInfo.EvoWeaver}
\alias{GeneVector.EvoWeaver}
\alias{Ancestral.EvoWeaver}
\docType{data}
\title{
Sequence-Level Predictions for EvoWeaver
}
\description{
\code{EvoWeaver} incorporates four classes of prediction, each with multiple
methods and algorithms. Sequence-Level (SL) methods examine conservation of patterns in sequence data, commonly exhibited due to physical interactions
between proteins.

\code{predict.EvoWeaver} currently supports three SL methods:
\itemize{
  \item \code{'SequenceInfo'}
  \item \code{'GeneVector'}
  \item \code{'Ancestral'}
}
}

\details{
All residue methods require a \code{EvoWeaver} object initialized
with \code{dendrogram} objects and ancestral states. See \code{\link{EvoWeaver}} for more information on input data types.

When \code{Method='Ensemble'} or \code{Method="SequenceLevel"}, EvoWeaver uses
methods \code{SequenceInfo} and \code{GeneVector}.

The \code{SequenceInfo} method looks at mutual information between sites in a multiple sequence alignment (MSA). This approach extends prior work in Martin et al. (2005). Each site from the first gene group is paired with the site from the second gene group that maximizes their mutual information.

The \code{GeneVector} method uses the natural vector encoding method introduced in Zhao et al. (2022). This encodes each gene sequences as a 92-dimensional vector, with the following entries:

\deqn{N(S) = (n_A,n_C,n_G,n_T,\\
              \qquad\qquad\;\,\mu_A,\mu_C,\mu_G,\mu_T,\\
              \qquad\qquad\quad\, D_2^A,D_2^C,D_2^G,D_2^T,\\
              \qquad\qquad\qquad n_{AA},n_{AC},\dots,n_{TT},\\
              \qquad\qquad\qquad\quad\;\; n_{AAA},n_{AAC},\dots,n_{TTT})}

Here \eqn{n_X} is the raw total count of nucleotide \eqn{X} (or di/trinucleotide). For single nucleotides, we also calculate \eqn{\mu_X}, the mean location of nucleotide \eqn{X}, and \eqn{D_2^X}, the second moment of the location of nucleotide \eqn{X}. The overall natural vector for a COG is calculated as the normalized mean vector from the natural vectors of all component gene sequences. Interaction scores are computed using Pearson's R between each COG's natural vector. These di/trinucleotide counts are by default excluded, but can be included using the \code{extended=TRUE} argument. Using the extended counts has shown minimal increased accuracy at the cost of slower runtime in benchmarking.

The \code{Ancestral} method calculates coevolution by looking at correlation of residue mutations near the leaves of each respective gene tree.
}

\format{
  None.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\references{
Martin, L. C., Gloor, G. B., Dunn, S. D. & Wahl, L. M, \emph{Using information theory to search for co-evolving residues in proteins}. Bioinformatics, 2005. \bold{21}(4116-4124).

Zhao, N., et al., \emph{Protein-protein interaction and non-interaction predictions using gene sequence natural vector}. Nature Communications Biology, 2022. \bold{5}(652).
}

\seealso{

\code{\link{EvoWeaver}}

\code{\link{predict.EvoWeaver}}

\link[=EvoWeaver-PPPreds]{EvoWeaver Phylogenetic Profiling Predictors}

\link[=EvoWeaver-PSPreds]{EvoWeaver Phylogenetic Structure Predictors}

\link[=EvoWeaver-GOPreds]{EvoWeaver Gene Organization Predictors}
}
