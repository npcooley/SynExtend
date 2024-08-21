\name{DecisionTree-class}
\alias{DecisionTree-class}
\alias{plot.DecisionTree}
\alias{as.dendrogram.DecisionTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Decision Trees for Random Forests
}
\description{
DecisionTree objects comprising random forest models generated with \code{\link{RandForest}}.
}
\usage{
\method{as.dendrogram}{DecisionTree}(object, ...)

\method{plot}{DecisionTree}(x, plotGain=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  an object of class \code{DecisionTree} to convert to class \code{\link{dendrogram}}.
}
  \item{x}{
  an object of class \code{DecisionTree} to plot.
  }
  \item{plotGain}{
  logical; should the Gini gain or decrease in sum of squared error be plotted for each decision point of the tree?
  }
  \item{...}{
  For \code{plot}, further arguments passed to \code{\link{plot.dendrogram}} and \code{\link{text}}. Arguments prefixed with \code{"text."} (e.g., \code{text.cex}) will be passed to \code{text}, and all other arguments are passed to \code{plot.dendrogram}.

  For \code{as.dendrogram}, \code{...} is further arguments for consistency with the generic definition.
  }
}
\details{
These methods help work with \code{DecisionTree} objects, which are returned as part of \code{\link{RandForest}}. Coercion to \code{\link{dendrogram}} objects creates a \code{'dendrogram'} corresponding to the structure of the decision tree. Each internal node possesses the standard attributes present in a \code{'dendrogram'} object, along with the following extra attributes:

\itemize{
\item \code{variable}: which variable was used to split at this node.
\item \code{thresh}: cutoff for partitioning points; values less than \code{thresh} are assigned to the left node, and those greater than to the right node.
\item \code{gain}: change in the metric to maximize. For classification trees this is the Gini Gain, and for regression trees this is the decrease in sum of squared error.
}

Plotting allows for extra arguments to be passed to \code{plot} and \code{text}. Arguments prefixed with \code{'text'} are passed to \code{text}, which controls the labeling of internal nodes. Common arguments used here are \code{text.cex}, \code{text.adj}, \code{text.srt}, and \code{text.col}. All other arguments are passed to \code{plot.dendrogram}. For example, \code{col='blue'} would change the dendrogram color to blue, whereas \code{text.col='blue'} would change the interior node labels to blue (but not the dendrogram itself).
}
\value{
\code{as.dendrogram} returns an object of class \code{'dendrogram'}. \code{plot} returns \code{NULL} invisibly.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}
\section{Warning}{
These functions can be quite slow for large decision trees. Usage is discouraged for trees with more than 100 internal nodes.
}

\seealso{
\code{\link{RandForest}}
}
\examples{
set.seed(199L)
n_samp <- 100L
AA <- rnorm(n_samp, mean=1, sd=5)
BB <- rnorm(n_samp, mean=2, sd=3)
CC <- rgamma(n_samp, shape=1, rate=2)
err <- rnorm(n_samp, sd=0.5)
y <- AA + BB + 2*CC + err

d <- data.frame(AA,BB,CC,y)
train_i <- 1:90
test_i <- 91:100
train_data <- d[train_i,]
test_data <- d[test_i,]

rf_regr <- RandForest(y~., data=train_data, rf.mode="regression", max_depth=5L)
if(interactive()){
  # Visualize one of the decision trees
  plot(rf_regr[[1]])
}

dend <- as.dendrogram(rf_regr[[1]])
plot(dend)
}
