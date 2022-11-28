\name{ProtWeaver-PAPreds}
\alias{ProtWeaver-PAPreds}
\alias{Jaccard.ProtWeaver}
\alias{Hamming.ProtWeaver}
\alias{MutualInformation.ProtWeaver}
\alias{ProfDCA.ProtWeaver}
\alias{Behdenna.ProtWeaver}
\alias{GainLoss.ProtWeaver}
\alias{CorrGL.ProtWeaver}
\docType{data}
\title{
Presence/Absence Predictions for ProtWeaver
}
\description{
\code{ProtWeaver} incorporates four classes of prediction, each with multiple 
methods and algorithms. Presence/Absence (PA) methods examine conservation 
of gain/loss events within orthology groups using phylogenetic profiles
constructed from presence/absence patterns.

\code{predict.ProtWeaver} currently supports six PA methods:
\itemize{
  \item \code{'Jaccard'}
  \item \code{'Hamming'}
  \item \code{'CorrGL'}
  \item \code{'MutualInformation'}
  \item \code{'ProfDCA'}
  \item \code{'Behdenna'}
  \item \code{'GainLoss'}
}
}

\details{
Most PA methods are compatible with a \code{ProtWeaver} object initialized
with any input type. See \code{\link{ProtWeaver}} for more information on input data types.

All of these methods use PA profiles, which are binary presence/absence
vectors such that 1 implies the corresponding genome has that particular
gene, and 0 implies the genome does not have that particular gene.

Methods \code{Hamming} and \code{Jaccard} use Hamming and Jaccard distance
(respectively) of PA profiles to determine overall score.

\code{MutualInformation} uses mutual information of PA profiels to determine
score, employing a weighting scheme such that \eqn{11} and \eqn{00} give
positive information, and \eqn{10} and \eqn{01} give negative information. 

\code{ProfDCA} uses the direct coupling analysis algorithm introduced by
Weigt et al. (2005) to determine direct information between PA profiles.
This approach has been validated on PA profiles in Fukunaga and Iwasaki (2022),
though the implementation in \code{ProtWeaver} forsakes the persistent contrasive divergence method in favor of the the algorithm from
Lokhov et al. (2018) for increased speed and exact solutions.

\code{Behdenna} implements the method detailed in Behdenna et al. (2016) to
find statistically significant interactions using co-occurence of gain/loss 
events mapped to ancestral states on a species tree. This method requires
a species tree as input. If the \code{ProtWeaver} object is initialized with \code{dendrogram} objects, \code{\link{SuperTree}} will be used to infer a species tree.

\code{GainLoss} uses a similar method to \code{Behdenna}. This method uses Fitch Parsimony to infer where events were gained or lost on a species tree, and then looks for distance between these gain/loss events. Unlike \code{Behdenna}, this method takes into account the types of events (ex. gain/gain and loss/loss are treated differently than gain/loss). This method requires
a species tree as input. If the \code{ProtWeaver} object is initialized with \code{dendrogram} objects, \code{\link{SuperTree}} will be used to infer a species tree.

\code{CorrGL} infers where events were gained or lost on a species tree as in method
\code{GainLoss}, then uses a Pearson's correlation coefficient weighted by p-value to infer similarity.
}

\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\references{
Behdenna, A., et al., \emph{Testing for Independence between Evolutionary Processes.} Systematic Biology, 2016. \bold{65}(5): p. 812-823.

Date, S.V. and E.M. Marcotte, \emph{Discovery of uncharacterized cellular systems by genome-wide analysis of functional linkages.} Nature Biotechnology, 2003. \bold{21}(9): p. 1055-1062.

Fukunaga, T. and W. Iwasaki, \emph{Inverse Potts model improves accuracy of phylogenetic profiling.} Bioinformatics, 2022.

Lokhov, A.Y., et al., \emph{Optimal structure and parameter learning of Ising models.} Science advances, 2018. \bold{4}(3): p. e1700791.

Pellegrini, M., et al., \emph{Assigning protein function by comparative genome analysis: Protein phylogenetic profiles.} Proceedings of the National Academy of Sciences, 1999. \bold{96}(8) p. 4285-4288

Weigt, M., et al., \emph{Identification of direct residue contacts in protein-protein interaction by message passing.} Proceedings of the National Academy of Sciences, 2009. \bold{106}(1): p. 67-72.
}

\seealso{

\code{\link{ProtWeaver}}

\code{\link{predict.ProtWeaver}}

\link[=ProtWeaver-DMPreds]{ProtWeaver Distance Matrix Predictors}

\link[=ProtWeaver-ColocPreds]{ProtWeaver Co-localization Predictors}

\link[=ProtWeaver-ResiduePreds]{ProtWeaver Residue Level Predictors}
}