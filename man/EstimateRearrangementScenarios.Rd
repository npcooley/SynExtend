\name{EstimateRearrangementScenarios}
\alias{EstimateRearrangementScenarios}

\title{
Estimate Genome Rearrangement Events with Double Cut and Join Operations
}
\description{
Take in a Synteny object from FindSynteny and return predicted rearrangement events.
}
\usage{
EstimateRearrangementScenarios(synt,
                               num_runs = -1, 
                               verbose = TRUE,
                               mean = FALSE,
                               opp_value = 0,
                               min_unique_length = 10,
                               min_block_length = -1,
                               actual = c(),
                               test_run = 0)
}

\arguments{
  \item{synt}{
  Synteny object obtained from running the FindSynteny object. Expected input is 
  unichromosomal sequences, though multichromosomal sequences are supported.
  }
  \item{num_runs}{
  Numeric; Number of times to simulate scenarios. The default value of -1 (and all non-positive
  values) runs each analysis for \emph{2b} iterations, where \emph{b} is the number of unique
  breakpoints.
  }
  \item{mean}{
  If TRUE, returns the mean number of inversions and transpositions found. If FALSE,
  returns the scenario corresponding to the minimum total number of operations 
  across all runs. This parameter only affects the number of inversions and 
  transpositions reported; the specific scenario returned is one of the runs that 
  resulted in a minimum value.
  }
  \item{min_block_length}{
  Minimum size of syntenic blocks to use for analysis. The default value accepts all blocks.
  Set to a larger value to ignore sections of short mutations that could be the 
  result of SNPs or other small-scale mutations.
  }
  \item{verbose}{
  Displays progress bars if verbose=TRUE.
  }
  \item{opp_value}{
  A parameter.
  }
  \item{min_unique_length}{
  A parameter.
  }
  \item{actual}{
  A parameter.
  }
  \item{test_run}{
  A parameter.
  }
}
\details{
\code{EstimateRearrangementScenarios} is an implementation of the Double Cut and Join
(DCJ) method for analyzing large scale mutation events. This implementation incentivizes
minimum number of total events rather than total number of DCJs. Multiple scenarios are
computed, and results of estimation are returned to the user. 
}
\value{
An \emph{NxN} matrix with the same shape as the input Synteny object. 

The diagonal corresponds to total sequence length of the corresponding genome. 

In the upper triangle, entry
\code{[i,j]} corresponds to the percent hits between genome \code{i} and genome \code{j}.
In the lower triangle, entry \code{[i,j]} contains a List object with 5 properties:
\itemize{
\item \code{$Inversions} and \code{$Transpositions} contain the (mean/min) number of estimated inversions
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
}
\examples{
db <- system.file("extdata", "Influenza.sqlite", package="DECIPHER")
synteny <- FindSynteny(db)
synteny

rearrs <- EstimateRearrangementScenarios(synteny)

rearrs          # view whole object
rearrs[[2,1]]   # view details on Genomes 1 and 2
}

