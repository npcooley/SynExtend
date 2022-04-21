\name{EstimRearrScen}
\alias{EstimRearrScen}
\alias{EstimateRearrangementScenarios}

\title{
Estimate Genome Rearrangement Events with Double Cut and Join Operations
}
\description{
Take in a \code{\link[=Synteny-class]{Synteny}} object and return predicted rearrangement events.
}
\usage{
EstimRearrScen(SyntenyObject, NumRuns = -1,
                Mean = FALSE, MinBlockLength = -1,
                Verbose = TRUE)
}

\arguments{
  \item{SyntenyObject}{
  \code{\link[=Synteny-class]{Synteny}} object, as obtained from running \code{\link{FindSynteny}}. 
  Expected input is unichromosomal sequences, though multichromosomal sequences are supported.
  }
  \item{NumRuns}{
  Numeric; Number of times to simulate scenarios. The default value of -1 (and all non-positive
  values) runs each analysis for \eqn{\sqrt b} iterations, where \code{b} is the number of unique
  breakpoints.
  }
  \item{Mean}{
  Logical; If TRUE, returns the mean number of inversions and transpositions found. If FALSE,
  returns the scenario corresponding to the minimum total number of operations 
  across all runs. This parameter only affects the number of inversions and 
  transpositions reported; the specific scenario returned is one of the runs that 
  resulted in a minimum value.
  }
  \item{MinBlockLength}{
  Numeric; Minimum size of syntenic blocks to use for analysis. The default value accepts all blocks.
  Set to a larger value to ignore sections of short mutations that could be the 
  result of SNPs or other small-scale mutations.
  }
  \item{Verbose}{
  Logical; indicates whether or not to display a progress bar and print the time difference upon completion.
  }
}
\details{
\code{EstimRearrScen} is an implementation of the Double Cut and Join
(DCJ) method for analyzing large scale mutation events. 

The DCJ model is commonly used to model genome rearrangement operations. Given a genome,
we can create a connected graph encoding the order of conserved genomic regions. 
Each syntenic region is split into two nodes, with one encoding the beginning and one encoding the end
(beginning and end defined relative to the direction of transcription). Each node
is then connected to the two nodes it is adjacent to in the genome.

For example, given a genome with 3 syntenic regions \eqn{a-b-c} such that \eqn{b}
is transcribed in the opposite direction relative to \eqn{a,c}, our graph would consist 
of nodes and edges \eqn{a1-a2-b2-b1-c1-c2}.

Given two genomes, we derive syntenic regions between the two samples and then construct
two of these graph structures. A DCJ operation is one that cuts two connections 
of a common color and creates two new edges. The goal of the DCJ model is to rearrange
the graph of the first genome into the second genome using DCJ operations. The
DCJ distance is defined as the minimum number of DCJ operations to transform one
graph into another.

It can be easily shown that inversions can be performed with a single DCJ operation,
and block interchanges/order rearrangements can be performed with a sequence of 
two DCJ operations. DCJ distance defines a metric space, and prior work has 
demonstrated algorithms for fast computation of the DCJ distance. 

However, DCJ distance inherently incentivizes inversions over block interchanges 
due to the former requiring half as many DCJ operations. This is a strong assumption,
and there is no evidence to support gene order rearrangements occuring half as often
as gene inversions.

This implementation incentivizes minimum number of total events rather than total number of DCJs. 
As the search space is large and multiple sequences of events can be equally parsimonious,
this algorithm computes multiple scenarios with random sequences of operations to
try to find the minimum amount of events. Users can choose to receive the best found
solution or the mean number of events from all solutions.
}
\value{
An \emph{NxN} matrix with the same shape as the input Synteny object. 

The diagonal corresponds to total sequence length of the corresponding genome. 

In the upper triangle, entry
\code{[i,j]} corresponds to the percent hits between genome \code{i} and genome \code{j}.
In the lower triangle, entry \code{[i,j]} contains a List object with 5 properties:
\itemize{
\item \code{$Inversions} and \code{$Transpositions} contain the (Mean/min) number of estimated inversions
and transpositions (resp.) between genome \code{i} and genome \code{j}. 
\item \code{$pct_hits}
contains percent hits between the genomes. 
\item \code{$Scenario} shows a rearrangement
scenario resulting in a minimum number of events between the genomes, using 
permutation matrix notation. Each line shows the operation performed, the permutation
vector corresponding to the current state relative to genome 1, and the number of 
blocks modified by the previous step.
\item \code{$Key} displays a key for converting from permutation matrix to genomes.
This matrix contains the start index for each block on each genome, the length of
the block, the relative direction of the block on genome 2 compared to genome 1, 
and the permutation number for that block.
}
}
\references{
Friedberg, R., Darling, A. E., & Yancopoulos, S. (2008). Genome rearrangement by the double cut and join operation. \emph{Bioinformatics}, 385-416.
}
\author{
Aidan Lakshman (\email{ahl27@pitt.edu})
}

\seealso{
\code{\link[DECIPHER]{FindSynteny}}

\code{\link[=Synteny-class]{Synteny}}
}
\examples{
db <- system.file("extdata", "Influenza.sqlite", package="DECIPHER")
synteny <- FindSynteny(db)
synteny

rearrs <- EstimRearrScen(synteny)

rearrs          # view whole object
rearrs[[2,1]]   # view details on Genomes 1 and 2
}

