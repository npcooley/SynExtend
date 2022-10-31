\name{GeneralizedRF}
\alias{GeneralizedRF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Information-Theoretic Generalized Robinson-Foulds Similarity
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GeneralizedRF(dend1, dend2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dend1}{
  An object of class \code{dendrogram}, representing a bifurcating phylogenetic
  tree.
}
  \item{dend2}{
  An object of class \code{dendrogram}, representing a bifurcating phylogenetic 
  tree.
}
}
\details{
This function implements Robinson-Foulds distance using the clustering information
approach first described in Smith (2020). This function iteratively pairs
internal tree branches of a phylogeny based on their similarity, then scores 
overall similarity as the sum of these measures. The similarity score is then
converted to a distance by normalizing by the average entropy of the two trees.
This metric has been demonstrated to outperform numerous other metrics in capabilities; 
see the original publication cited in References for more information. 

This implementation uses a quick pairing algorithm to speed up computation at 
the cost of slight inaccuracy. Results should be suitably close to the optimal
pairing in the vast majority of cases.

Note that this function requires the input dendrograms to be labeled alike (ex.
leaf labeled \code{abc} in \code{dend1} represents the same species as 
leaf labeled \code{abc} in \code{dend2}). 
Labels can easily be modified using \code{\link{dendrapply}}. 
}
\value{
A value from 0 to 1, with 1 indicating the trees are maximally different
and 0 indicating the trees are identical.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dend1, dend2) 
{
    stopifnot(`inputs must both be dendrograms!` = is(dend1, 
        "dendrogram") && is(dend2, "dendrogram"))
    incommonLabs <- intersect(labels(dend1), labels(dend2))
    if (length(incommonLabs) == 0) 
        return(0)
    tree1ptr <- .Call("initCDend", dend1)
    on.exit(rm(tree1ptr))
    tree2ptr <- .Call("initCDend", dend2)
    on.exit(rm(tree2ptr))
    val <- .Call("GRFInfo", incommonLabs)
    return(val)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
