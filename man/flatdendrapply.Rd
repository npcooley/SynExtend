\name{flatdendrapply}
\alias{flatdendrapply}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\code{sapply} Equivalent for Dendrogram Objects
}
\description{
This function allows for applying a function to dendrograms recursively, 
returning output simplified to vector if possible. This also allows for 
applying different functions to nodes and leaves, and allows for selectively 
applying a function to only internal nodes or leaves.
}
\usage{
flatdendrapply(dend, 
                NODEFUN=NULL, 
                LEAFFUN=NODEFUN, 
                INCLUDEROOT=TRUE, ...)
}
\arguments{
  \item{dend}{ 
  an object of class \code{\link{dendrogram}}
  }
  
  \item{NODEFUN}{
  A function to apply to each internal node. If \code{NULL}, internal nodes are 
  skipped.
  }
  
  \item{LEAFFUN}{
  A function to apply to each leaf node. If \code{NULL}, leaf nodes are skipped.
  Defaults to \code{NODEFUN}, meaning \code{NODEFUN} will be applied to all nodes 
  (interior and leaves) if \code{LEAFFUN} is unspecified.
  }
  
  \item{INCLUDEROOT}{
  If \code{FALSE}, skips applying \code{NODEFUN} to the root node.
  }
  
  \item{...}{
  Additional arguments passed to BOTH \code{NODEFUN} and \code{LEAFFUN}
  }
}
\details{
\code{\link{flatdendrapply}} applies \code{NODEFUN} to all interior nodes
and \code{LEAFFUN} to all leaf nodes, then returns the result as a flat list
(a list with no nested components). If the result can be simplified to a vector
without loss of information (meaning all components are atomic and length 1), 
the function will return a vector. 

This function contrasts with \code{\link{dendrapply}}, which applies a single
function to all nodes of the dendrogram and returns a nested list with identical
topology to the input dendrogram. \code{flatdendrapply} sacrifices the original
dendrogram structure to return input that is easier to work with (since return
values are guaranteed to be a list with no nested elements). \code{flatdendrapply}
is used internally in some other functions.
}
\value{
Either a list or a vector. 

If the function applied to each node/leaf returns an
atomic value with length 1, then the return value is a vector of length \code{n},
with \code{n} the number of nodes (or leaves/interior nodes, if \code{NODEFUN}
or \code{LEAFFUN} are \code{NULL}, respectively).

If the function applied to each node/leaf returns anything other than an atomic 
value with length 1, then the return value is a list of with the same length as 
described above. This list is guaranteed to not have any nested lists.
}

\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\seealso{
\code{\link{dendrogram}}, \code{\link{dendrapply}}, \code{\link{rapply}}
}
\examples{
data("SuperTreeEx")

exDend <- exData[[1]]

# Get labels of all leaf nodes of a dendrogram
flatdendrapply(exDend, LEAFFUN = attr, which='label')

# Get heights of all interior nodes
flatdendrapply(exDend, NODEFUN = attr, which='height', LEAFFUN = NULL)

# Get heights of all nodes (LEAFFUN = NODEFUN by default)
flatdendrapply(exDend, NODEFUN = attr, which='height')

# Get heights of interior nodes and labels of leaf nodes
flatdendrapply(exDend, NODEFUN = function(x) attr(x, 'height'),
                       LEAFFUN = function(x) attr(x, 'label'))
}
