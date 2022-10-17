\name{SuperTree}
\alias{SuperTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a Species Tree from Gene Trees
}
\description{
Given a set of unrooted gene trees, creates a species tree. This function
works for rooted gene trees, but may not accurately root the resulting tree.
}
\usage{
SuperTree(myDendList, NAMEFUN=NULL, Verbose=TRUE, Processors=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{myDendList}{
  List of \code{dendrogram} objects, where each entry is an unrooted gene tree.
}
\item{NAMEFUN}{
  Optional input specifying a function to apply to each leaf to convert gene tree leaf 
  labels into species names. This function should take as input a character vector
  and return a character vector of the same size. By default equals \code{NULL}, 
  indicating that gene tree leaves are already labeled with species identifiers.
  See details for more information.
}

\item{Verbose}{
Should output be displayed?
}

\item{Processors}{
Number of processors to use for calculating the final species tree.
}
}
\details{
This implementation follows the ASTRID algorithm for estimating a species tree
from a set of unrooted gene trees. Input gene trees are not required to have
identical species sets, as the algorithm can handle missing entries in gene trees.
See the original paper linked in the references section for more information.

\code{SuperTree} allows an optional argument called \code{NAMEFUN} to apply a 
renaming step to leaf labels. Gene trees as constructed by other functions in 
\code{SynExtend} (ex. \code{\link{DisjointSet}}) often include other information 
aside from species name when labeling genes, but \code{SuperTree} requires that 
leaf nodes of the gene tree are labeled with just an identifier corresponding to
which species/genome each leaf is from. Duplicate values are allowed. See the examples
section for more details on what this looks like and how to handle it.
}
\value{
A \code{\link{dendrogram}} object corresponding to the species tree constructed
from input gene trees.
}
\references{
Vachaspati, P., Warnow, T. \emph{ASTRID: Accurate Species TRees from Internode Distances}. BMC Genomics, 2015. \bold{16} (Suppl 10): S3.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\seealso{
\code{\link[DECIPHER]{TreeLine}}, \code{\link{SuperTreeEx}}
}
\examples{
# Loads a list of dendrograms
# each is a gene tree from Streptomyces genomes
data("SuperTreeEx", package="SynExtend")

# Notice that the labels of the tree are in #_#_# format
# See the man page for SuperTreeEx for more info
labs <- labels(exData[[1]])
labs

# The first number corresponds to the species,
# so we need to trim the rest in each leaf label
namefun <- function(x) gsub("([0-9A-Za-z]*)_.*", "\\\\1", x)
namefun(labs) # trims to just first number

# This function replaces gene identifiers with species identifiers
# we pass it to NAMEFUN
# Note NAMEFUN should take in a character vector and return a character vector
tree <- SuperTree(exData, NAMEFUN=namefun)
}
