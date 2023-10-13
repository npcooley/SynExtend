\name{EvoWeaver-GOPreds}
\alias{EvoWeaver-GOPreds}
\alias{Coloc.EvoWeaver}
\alias{ColocMoran.EvoWeaver}
\alias{TranscripMI.EvoWeaver}
\docType{data}
\title{
Gene Organization Predictions for EvoWeaver
}
\description{
\code{EvoWeaver} incorporates four classes of prediction, each with multiple 
methods and algorithms. Co-localization (Coloc) methods examine
conservation of relative location and transcriptional direction of genetic
regions within the genome.

\code{predict.EvoWeaver} currently supports three Coloc methods:
\itemize{
  \item \code{'Coloc'}
  \item \code{'ColocMoran'}
  \item \code{'TranscripMI'}
}
}

\details{
All distance matrix methods require a \code{EvoWeaver} object initialized
with gene locations using the a three or four number code. See \code{\link{EvoWeaver}} for more information on input data types.

The built-in \code{Coloc} examines relative location of genes within genomes
as evidence of interaction. For a given pair of genes, the score is given by
\eqn{\sum_{G}e^{1-|dI_G|}}, where \eqn{G} the set of genomes and \eqn{dI_G} the 
difference in index between the two genes in genome \eqn{G}. Using gene index
instead of number of base pairs avoids bias introduced by gene and genome length.

\code{ColocMoran} measures the extent to which gene distances are preserved across a phylogeny. This function uses the same initial scoring scheme as \code{Coloc}, but can handle paralogs. The raw scores are passed into \code{\link{MoransI}} to calculate spatial autocorrelation. "Space" is taken as \eqn{e^{-C}}, where \eqn{C} is the Cophenetic distance matrix calculated from the species tree of the inputs. As such, this method requires a species tree as input, which can be calculated from a set of gene trees using \code{\link{SuperTree}}.

\code{TranscripMI} uses mutual information of the transcriptional direction of each pair of genes. Conservation of relative transcriptional direction between gene pairs has been shown to imply functional association in prior work. This algorithm requires that the \code{EvoWeaver} object is initialized with a four number code, with the third number either \code{0} or \code{1}, denoting whether the gene is on the forward or reverse strand. The mutual information is calculated as: 

\deqn{\sum_{x \in X}\sum_{y \in Y}(-1)^{(x!=y)}P_{(X,Y)}(x,y)\; \log\left(\frac{P_{(X,Y)}(x,y)}{P_X(x)P_Y(y)}\right)}

Here \eqn{X=Y=\{0,1\}}, \eqn{x} is the direction of the gene with lower index, \eqn{y} is the direction of the gene with higher index, and \eqn{P_{(T)}(t)} is the probability of \eqn{T=t}. Note that this is a weighted MI as introduced by Beckley and Wright (2021). The mutual information is augmented by the addition of a single pseudocount to each value, and normalized by the joint entropy of \eqn{X,Y}. P-values are calculated using Fisher's Exact Test on the contingency table. 
}

\value{
  None.
}

\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\references{
Beckley, Andrew and E. S. Wright. \emph{Identification of antibiotic pairs that evade concurrent resistance via a retrospective analysis of antimicrobial susceptibility test results}. The Lancet Microbe, 2021. \bold{2}(10): 545-554.

Korbel, J. O., et al., \emph{Analysis of genomic context: prediction of functional associations from conserved bidirectionally transcribed gene pairs}. Nature Biotechnology, 2004. \bold{22}(7): 911-917.

Moran, P. A. P., \emph{Notes on Continuous Stochastic Phenomena}. Biometrika, 1950. \bold{37}(1): 17-23.
}

\seealso{

\code{\link{EvoWeaver}}

\code{\link{predict.EvoWeaver}}

\link[=EvoWeaver-PPPreds]{EvoWeaver Phylogenetic Profiling Predictors}

\link[=EvoWeaver-PSPreds]{EvoWeaver Phylogenetic Structure Predictors}

\link[=EvoWeaver-SLPreds]{EvoWeaver Sequence-Level Predictors}
}