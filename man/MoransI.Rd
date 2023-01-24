\name{MoransI}
\alias{MoransI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Moran's \emph{I} Spatial Autocorrelation Index
}
\description{
Calculates Moran's \emph{I} to measure spatial autocorrelation for a set of signals dispersed in space.
}
\usage{
MoransI(values, weights, alternative='two.sided')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{values}{
  Numeric vector containing signals for each point in space.
}
\item{weights}{
  Distances between each point in space. This should be a numeric object of class \code{\link{dist}} with \code{Size} attribute equivalent to the length of \code{values}. 
}

\item{alternative}{
  For hypothesis testing against the null of no spatial correlation, how should a p-value be calculated? Should be one of \code{c("two.sided", "less", "greater")}, or an unambiguous abbreviation.
}
}
\details{
  Moran's \emph{I} is a measure of how much the spatial arrangement of a set of datapoints correlates with the value of each datapoint. The index takes a value in the range \eqn{[-1,1]}, with values close to 1 indicating high correlation between location and value (points have increasingly similar values as they increase in proximity), values close to -1 indicating anticorrelation(points have increasingly different values as they increase in proximity), and values close to 0 indicating no correlation.
  
  The value itself is calculated as:

\deqn{I = \frac{N}{W}\frac{\sum_i^N \sum_j^N w_{ij}(x_i - \bar x)(x_j - \bar x)}{\sum_i^N (x_i - \bar x)^2}}

Here, \eqn{N} is the number of points, \eqn{w_{ij}} is the distance between points \eqn{i} and \eqn{j}, \eqn{W = \sum_{i,j} w_{ij}} (the sum of all the weights), \eqn{x_i} is the value of point \eqn{i}, and \eqn{\bar x} is the sample mean of the values.

Moran's \emph{I} has a closed form calculation for variance and expected value, which are calcalated within this function. The full form of the variance is fairly complex, but all the equations are available for reference \href{https://en.wikipedia.org/wiki/Moran\%27s_I#Expected_value}{here}.

A p-value is estimated using the expected value and variance using a null hypothesis of no spatial autocorrelation, and the alternative hypothesis specified in the \code{alternative} argument. Note that if fewer than four datapoints are supplied, the variance of Moran's I is infinite. The function will return a standard deviation of \code{Inf} and a p-value of 1 in this case.
}
\value{
A \code{\link{list}} object containing the following named values:
  \itemize{
    \item{
     \code{observed}: The value of Moran's \emph{I} (numeric in the range \eqn{[-1,1]}).
    }
    \item{
     \code{expected}: The expected value of Moran's \emph{I} for the input data.
    }
    \item{
     \code{sd}: The standard deviation of Moran's \emph{I} for the input data.
    }
    \item{
     \code{p.value}: The p-value for the input data, calculated with the alternative hypothesis as specified in \code{alternative}.
    }
  }
}
\references{
Moran, P. A. P. \emph{Notes on Continuous Stochastic Phenomena}. Biometrika, 1950. \bold{37}(1): 17-23.

Gittleman, J. L. and M. Kot. \emph{Adaptation: Statistics and a Null Model for Estimating Phylogenetic Effects}. Systematic Zoology, 1990. \bold{39}:227-241.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\examples{
# Make a distance matrix for a set of 50 points
# These are just random numbers in the range [0.1,2]
NUM_POINTS <- 50
distmat <- as.dist(matrix(runif(NUM_POINTS**2, 0.1, 2), 
                          ncol=NUM_POINTS))

# Generate some random values for each of the points
vals <- runif(NUM_POINTS, 0, 3) 

# Calculate Moran's I
MoransI(vals, distmat, alternative='two.sided')

# effect size should be pretty small
# and p-value close to 0.5
# since this is basically random data
}
