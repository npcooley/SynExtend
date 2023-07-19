\name{FitchParsimony}
\alias{FitchParsimony}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate ancestral states using Fitch Parsimony
}
\description{
Ancestral states for binary traits can be inferred from presence/absence patterns at the tips of a dendrogram using Fitch Parsimony. This function works for an arbitrary number of states on bifurcating \code{dendrogram} objects.
}
\usage{
FitchParsimony(dend, num_traits, traits_list,
                  initial_state=rep(0L,num_traits), 
                  fill_ambiguous=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dend}{
  An object of class \code{'dendrogram'}
}
  \item{num_traits}{
  The number of traits to inferred, as an integer.
}
  \item{traits_list}{
  A list of character vectors, where the \code{i}'th entry corresponds to the leaf labels that have the trait \code{i}. 
}
  \item{initial_state}{
  The state assumed for the root node. Set to \code{NULL} to disable autofilling the root state.
  }
  \item{fill_ambiguous}{
  If \code{TRUE}, states that remain ambiguous after completion of the algorithm are filled in randomly.
  }
}
\details{
Fitch Parismony allows for fast inference of ancestral states of binary traits. The algorithm proceeds in three steps. 

First, traits are inferred upwards based on child nodes. If the child nodes have the same state (\code{1/1} or \code{0/0}), then the parent node is also set to that state. If the states are different, the parent node is set to \code{2}, denoting an ambiguous entry. If one child is ambiguous and the other is not, the parent is set to the non-ambiguous entry.

Second, traits are inferred downward to attempt to fill in ambiguous entries. If a node is not ambiguous but its child is, the child's state is set to the parent state. If specified, the root node's state is set to \code{initial_state} prior to this step.

Third, traits that remain ambiguous are optionally filled in (only if \code{fill_ambiguous} is set to \code{TRUE}). This proceeds by randomly setting ambiguous traits to either \code{1} or \code{0}.

The result is stored in the \code{FitchState} attribute within each node. 
}
\value{
A dendrogram with attribute \code{FitchState} set for each node, where this attribute is a binary vector of length \code{num_traits}.
}
\note{
It's FitchParsimony because this implementation is entirely in R, as opposed to internal \code{SynExtend} methods that utilize a slightly faster C-based implementation that is not user-exposed.
}
\references{
Fitch, Walter M. \emph{Toward defining the course of evolution: minimum change for a specific tree topology.} Systematic Biology, 1971. \bold{20}(4): p. 406-416.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\examples{
d <- as.dendrogram(hclust(dist(USArrests), "ave"))
labs <- labels(d)

# Defining some presence absence patterns
set.seed(123L)
pa_1 <- sample(labs, 15L)
pa_2 <- sample(labs, 20L)

# inferring ancestral states
fpd <- FitchParsimony(d, 2L, list(pa_1, pa_2))

# Checking a state
attr(fpd[[1L]], 'FitchState')

# Visualizing the results for the first pattern
# Tips show P/A patterns, edges show gain/loss (green/red)
fpd <- dendrapply(fpd, \(x){
  ai <- 1L
  s <- attr(x, 'FitchState')
  l <- list()
  
  if(is.leaf(x)){
    # coloring tips based presence/absence
    l$col <- ifelse(s[ai]==1L, 'green', 'red')
    l$pch <- 19
    attr(x, 'nodePar') <- l
  } else {
    # coloring edges based on gain/loss
    for(i in seq_along(x)){
      sc <- attr(x[[i]], 'FitchState')
      if(s[ai] != sc[ai]){
        l$col <- ifelse(s[ai] == 1L, 'red', 'green')
      } else {
        l$col <- 'black'
      }
      attr(x[[i]], 'edgePar') <- l
    }
  }
  
  x
}, how='post.order')
plot(fpd, leaflab='none')
}

