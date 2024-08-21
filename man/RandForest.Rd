\name{RandForest}
\alias{RandForest}
\alias{predict.RandForest}
\alias{RandForest.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Classification and Regression with Random Forests
}
\description{
\code{RandForest} implements a version of Breiman's random forest algorithm for classification and regression.
}
\usage{
RandForest(formula, data, subset, verbose=interactive(),
           weights, na.action,
           method='rf.fit',
           rf.mode=c('auto', 'classification', 'regression'),
           contrasts=NULL, ...)

\method{predict}{RandForest}(object, newdata=NULL,
                na.action=na.pass, ...)

## Called internally by `RandForest`
RandForest.fit(x, y=NULL,
    verbose=interactive(), ntree=10,
    mtry=floor(sqrt(ncol(x))),
    weights=NULL, replace=TRUE,
    sampsize=if(replace) nrow(x) else ceiling(0.632*nrow(x)),
    nodesize=1L, max_depth=NULL,
    method=NULL,
    terms=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    an object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. See \code{\link{lm}} for more details.
}
  \item{data}{
  An optional data frame, list, or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{RandForest} is called.
  }
  \item{subset}{
  an optional vector specifying a subset of observations to be used in the fitting process.
  }
  \item{weights}{
  an optional vector of weights to be used in the fitting process. Should be \code{NULL} or a numeric vector.
  }
  \item{na.action}{
  a function which indicates what should happen when the data contain \code{NAs}. Currently experimental.
  }
  \item{method}{
  currently unused.
  }
  \item{rf.mode}{
  one of \code{"auto"}, \code{"classification"}, \code{"regression"} (or an unambiguous abbreviation), specifying the type of trees to build. If \code{rf.mode="auto"}, the mode is inferred based on the type of the response variable.
  }
  \item{contrasts}{
  currently experimental; see \code{\link{lm}}.
  }
  \item{...}{
  further arguments passed to \code{RandForest.fit}.
  }
  \item{object}{
  an object of class \code{'RandForest'} for prediction.
  }
  \item{newdata}{
  new data to predict on, typically provided as a \code{data.frame} object.
  }
  \item{verbose}{
  logical: should progress be displayed?
  }
  \item{ntree}{
  number of decision trees to grow.
  }
  \item{mtry}{
  number of variables to try at each split.
  }
  \item{replace}{
  logical; should data be sampled with replacement during training?
  }
  \item{sampsize}{
  number of datapoints to sample for training each component decision tree.
  }
  \item{nodesize}{
  number of datapoints to stop classification (see "Details")
  }
  \item{max_depth}{
  maximum depth of component decision trees.
  }
  \item{x}{
    used internally by \code{RandForest.fit}
  }
  \item{y}{
    used internally by \code{RandForest.fit}
  }
  \item{terms}{
    used internally by \code{RandForest.fit}
  }
}
\details{
\code{RandForest} implements a version of Breiman's original algorithm to train a random forest model for classification or regression. Random forests are comprised of a set of decision trees, each of which is trained on a subset of the available data. These trees are individually worse predictors than a single decision tree trained on the entire dataset. However, averaging predictions across the ensemble of trees forms a model that is often more accurate than single decision trees while being less susceptible to overfitting.

Random forests can either be trained for classification or regression. Classification forests are comprised of trees that assign inputs to a specific class. The output prediction is a vector comprised of the proportion of trees in the forest that assigned the datapoint to each available class. Regresssion forests are comprised of trees that assign each datapoint to a single continuous value, and the output prediction is comprised of the mean prediction across all component trees. When \code{rf.mode="auto"}, the random forest will be trained in classification mode for response of type \code{"factor"}, and in regression mode for response of type \code{"numeric"}.

Several parameters exist to tune the behavior of random forests. The \code{ntree} argument controls how many decision trees are trained. At each decision point, the decision trees consider a random subset of available variables--the number of variables to sample is controlled by \code{mtry}. Each decision tree only sees a subset of available data to reduce its risk of overfitting. This subset is comprised of \code{sampsize} datapoints, which are sampled with or without replacement according to the \code{replace} argument.

Finally, the default behavior is to grow decision trees until they have fully classified all the data they see for training. However, this may lead to overfitting. Decision trees can be limited to smaller sizes by specifying the \code{max_depth} or \code{nodesize} arguments. \code{max_depth} refers to the depth of the decision tree. Setting this value to \code{n} means that every path from the root node to a leaf node will be at most length \code{n}. \code{nodesize} can be used to instead stop growing trees based on the size of the data to be partitioned at each decision tree node. If \code{nodesize=n}, then if a decision point receives less than \code{n} samples, it will stop trying to further split the data.

Classification forests are trained by maximizing the Gini Gain at each interior node. Split points are determined with exhaustive search for small data sizes, or simulated annealing for larger sizes. Regression forests are trained by maximizing the decrease in sum of squared error (SSE) if all points in each partition are assigned their mean output value. Nodes stop classification when either no partition improves the maximization metric (Gini Gain or decrease in SSE) or when the criteria specified by \code{nodesize} / \code{max_depth} are met.

Some of the arguments provided are for consistency with the base \code{\link{lm}} function. Use caution changing any values referred to as "Experimental" above. \code{NA} values may cause unintended behavior.
}
\value{
An object of class \code{'RandForest'}, which itself contains a number of objects of class \code{\link[=DecisionTree-class]{'DecisionTree'}} which can be used for prediction with \code{\link{predict.RandForest}}
}

\note{
Generating a single decision tree model is possible by setting \code{ntree=1} and \code{sampsize=nrow(data)}. \code{'DecisionTree'} objects do not currently support prediction.
}

\references{
Breiman, L. (2001), \emph{Random Forests}, Machine Learning 45(1), 5-32.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\seealso{
\link[=DecisionTree-class]{DecisionTree class}
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

## classification
y1 <- y < -5
y2 <- y < 0 & y >= -5
y3 <- y < 5 & y >= 0
y4 <- y >= 5
y_cl <- rep(0L, length(y))
y[y1] <- 1L
y[y2] <- 2L
y[y3] <- 3L
y[y4] <- 4L
d$y <- as.factor(y)
train_data <- d[train_i,]
test_data <- d[test_i,]

rf_classif <- RandForest(y~., data=train_data, rf.mode="classification", max_depth=5L)
if(interactive()){
  # Visualize one of the decision trees for classification
  plot(rf_classif[[1]])
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
