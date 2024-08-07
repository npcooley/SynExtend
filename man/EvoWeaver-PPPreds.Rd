\name{EvoWeaver-PPPreds}
\alias{EvoWeaver-PPPreds}
\alias{ExtantJaccard.EvoWeaver}
\alias{Hamming.EvoWeaver}
\alias{GLMI.EvoWeaver}
\alias{ProfDCA.EvoWeaver}
\alias{Behdenna.EvoWeaver}
\alias{GLDistance.EvoWeaver}
\alias{CorrGL.EvoWeaver}
\alias{PAJaccard.EvoWeaver}
\alias{PAOverlap.EvoWeaver}
\docType{data}
\title{
Phylogenetic Profiling Predictions for EvoWeaver
}
\description{
\code{EvoWeaver} incorporates four classes of prediction, each with multiple
methods and algorithms. Phylogenetic Profiling (PP) methods examine conservation
of gain/loss events within orthology groups using phylogenetic profiles
constructed from presence/absence patterns.

\code{predict.EvoWeaver} currently supports nine PP methods:
\itemize{
  \item \code{'ExtantJaccard'}
  \item \code{'Hamming'}
  \item \code{'GLMI'}
  \item \code{'PAPV'}
  \item \code{'CorrGL'}
  \item \code{'ProfDCA'}
  \item \code{'Behdenna'}
  \item \code{'GLDistance'}
  \item \code{'PAJaccard'}
  \item \code{'PAOverlap'}
}
}

\details{
Most PP methods are compatible with a \code{EvoWeaver} object initialized
with any input type. See \code{\link{EvoWeaver}} for more information on input data types.

When \code{Method='Ensemble'} or \code{Method="PhylogeneticProfiling"}, EvoWeaver uses
methods \code{GLMI}, \code{GLDistance}, \code{PAJaccard}, and \code{PAOverlap}.

All of these methods use presence/absence (PA) profiles, which are binary
vectors such that 1 implies the corresponding genome has that particular
gene, and 0 implies the genome does not have that particular gene.

Methods \code{Hamming} and \code{ExtantJaccard} use Hamming and Jaccard distance
(respectively) of PA profiles to determine overall score.

\code{GLMI} uses mutual information of gain/loss (G/L) vectors to determine
score, employing a weighting scheme such that concordant gains/losses give positive information,
discordant gains/losses give negative information, and events that do not cooccur with a gain/loss
in the other gene group give no information.

\code{PAJaccard} calculates the centered Jaccard index of P/A profiles, where each clade
with identical extant patterns is collapsed to a single leaf.

\code{PAOverlap} calculates the proportion of time in the ancestry that both genes cooccur
relative to the total time each individual gene occurs, based on ancestral states inferred
with Fitch parsimony.

\code{PAPV} calculates a p-value for PA profiles using Fisher's Exact Test. The returned score is provided as \code{1-p_value} so that larger scores indicate more significance, and smaller scores indicate less significance. This rescaling is consistent with the other similarity metrics in \code{EvoWeaver}. This can be used with \code{ExtantJaccard}, \code{Hamming}, or \code{GLMI} to weight raw scores by statistical significance.

\code{ProfDCA} uses the direct coupling analysis algorithm introduced by
Weigt et al. (2005) to determine direct information between PA profiles.
This approach has been validated on PA profiles in Fukunaga and Iwasaki (2022),
though the implementation in \code{EvoWeaver} forsakes the persistent contrasive divergence method in favor of the the algorithm from
Lokhov et al. (2018) for increased speed and exact solutions. Note that this algorithm is still extremely slow relative to the other methods despite the aforementioned runtime improvements.

\code{Behdenna} implements the method detailed in Behdenna et al. (2016) to
find statistically significant interactions using co-occurence of gain/loss
events mapped to ancestral states on a species tree. This method requires
a species tree as input. If the \code{EvoWeaver} object is initialized with \code{dendrogram} objects, \code{\link{SuperTree}} will be used to infer a species tree.

\code{GLDistance} uses a similar method to \code{Behdenna}. This method uses Fitch Parsimony to infer where events were gained or lost on a species tree, and then looks for distance between these gain/loss events. Unlike \code{Behdenna}, this method takes into account the types of events (ex. gain/gain and loss/loss are treated differently than gain/loss). This method requires
a species tree as input. If the \code{EvoWeaver} object is initialized with \code{dendrogram} objects, \code{\link{SuperTree}} will be used to infer a species tree.

\code{CorrGL} infers where events were gained or lost on a species tree as in method
\code{GLDistance}, then uses a Pearson's correlation coefficient weighted by p-value to infer similarity.
}

\value{
  None.
}

\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\references{
Behdenna, A., et al., \emph{Testing for Independence between Evolutionary Processes.} Systematic Biology, 2016. \bold{65}(5): p. 812-823.

Chung, N.C, et al., \emph{Jaccard/Tanimoto similarity test and estimation methods for biological presence-absence data.} BMC Bioinformatics, 2019. \bold{20}(S15).

Date, S.V. and E.M. Marcotte, \emph{Discovery of uncharacterized cellular systems by genome-wide analysis of functional linkages.} Nature Biotechnology, 2003. \bold{21}(9): p. 1055-1062.

Fukunaga, T. and W. Iwasaki, \emph{Inverse Potts model improves accuracy of phylogenetic profiling.} Bioinformatics, 2022.

Lokhov, A.Y., et al., \emph{Optimal structure and parameter learning of Ising models.} Science advances, 2018. \bold{4}(3): p. e1700791.

Pellegrini, M., et al., \emph{Assigning protein function by comparative genome analysis: Protein phylogenetic profiles.} Proceedings of the National Academy of Sciences, 1999. \bold{96}(8) p. 4285-4288

Weigt, M., et al., \emph{Identification of direct residue contacts in protein-protein interaction by message passing.} Proceedings of the National Academy of Sciences, 2009. \bold{106}(1): p. 67-72.
}

\seealso{

\code{\link{EvoWeaver}}

\code{\link{predict.EvoWeaver}}

\link[=EvoWeaver-PSPreds]{EvoWeaver Phylogenetic Structure Predictors}

\link[=EvoWeaver-GOPreds]{EvoWeaver Gene Organization Predictors}

\link[=EvoWeaver-SLPreds]{EvoWeaver Sequence-Level Predictors}
}