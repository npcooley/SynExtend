\name{subset.dendrogram}
\alias{subset.dendrogram}

\title{
Subsetting dendrogram objects
}
\description{
Subsets dendrogram objects based on leaf labels. Subsetting can either be by leaves to keep, or leaves to remove.

NOTE: This man page is specifically for \code{subset.dendogram}, see \code{?base::subset} for the generic \code{subset} function defined for vectors, matrices, and data frames.
}
\usage{
\method{subset}{dendrogram}(x, subset, invert=FALSE, ...)
}

\arguments{
  \item{x}{
  An object of class \code{'dendogram'}
}
\item{subset}{
  A vector of labels to keep (see \code{invert}).
}
\item{invert}{
  If set to \code{TRUE}, subset to only the leaves \emph{not} in \code{subset}.
}
\item{...}{
Additional arguments for consistency with generic.
}
}

\value{
An object of class \code{'dendrogram'} corresponding to the subsetted tree.
}

\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}
\note{
If none of the labels specified in the \code{subset} argument appear in the tree (or if all do when \code{invert=TRUE}), a warning is thrown and an empty object of class \code{'dendrogram'} is returned.
}

\seealso{
\code{\link[base]{subset}}
}
\examples{
d <- as.dendrogram(hclust(dist(USArrests), "ave"))

# Show original dendrogram
plot(d)

# Subset to first 10 labels
d1 <- subset(d, labels(d)[1:10])
plot(d1)

# Subset d1 to all except the first 2 labels
d2 <- subset(d1, labels(d1)[1:2], invert=TRUE)
plot(d2)
}

