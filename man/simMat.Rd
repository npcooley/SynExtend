\name{simMat}
\alias{simMat}
\alias{simMat-class}
\alias{as.simMat}
\alias{as.simMat.vector}
\alias{as.simMat.matrix}
\alias{Diag.simMat}
\alias{Diag<-.simMat}
\alias{Diag}
\alias{Diag<-}
\alias{print.simMat}
\alias{as.matrix.simMat}
\alias{as.data.frame.simMat}
\title{
Similarity Matrices
}
\description{
The \code{simMat} object is an internally utilized class that provides similar 
functionality to the \code{\link{dist}} object, but with matrix-like accessors.

Like \code{dist}, this object stores values as a vector, reducing memory by
making use of assumed symmetry. \code{simMat} currently only supports numeric
data types.
}
\usage{
## Create a blank sym object
simMat(VALUE, nelem, NAMES=NULL, DIAG=FALSE)

\method{as.simMat}{vector}(x, NAMES=NULL, DIAG=TRUE, ...)

\method{as.simMat}{matrix}(x, ...)

\method{print}{simMat}(x, n=10, ...)

\method{as.matrix}{simMat}(x, ...)

\method{as.data.frame}{simMat}(x, ...)

\method{Diag}{simMat}(x, ...)

\method{Diag}{simMat}(x) <- value
}
\arguments{
\item{VALUE}{
Numeric (or \code{NA_real_}) indicating placeholder values. A vector of values can
be provided for this function if desired.
}
\item{nelem}{
Integer; number of elements represented in the matrix. This corresponds to the
number of rows and columns of the object, so setting \code{nelem=10} would 
produce a \code{10x10} matrix.
}
\item{NAMES}{
Character (Optional); names for each row/column. If provided, this should be 
a character vector of length equal to \code{nelem}.
}
\item{DIAG}{
Logical; Is the diagonal included in the data? If \code{FALSE}, the constructor
generates 1s for the diagonal.
}
\item{x}{Various; for \code{print} and \code{Diag}, the \code{"simMat"} object to print. For \code{as.vector} or \code{as.matrix}, the vector or matrix (respectively). Note that \code{as.matrix} expects a symmetric matrix--providing a non-symmetric matrix will take only the upper triangle and produce a warning.}
\item{n}{Integer (Optional); maximum number of rows to print when printing a \code{simMat} object. Default is 10 rows.}
\item{value}{Numeric; value(s) to replace diagonal with.}
\item{...}{Additional parameters provided for consistency with generic.}
}
\details{
The \code{simMat} object has a very similar format to \code{dist} objects, but 
with a few notable changes:

\itemize{
\item \code{simMat} objects have streamlined \code{print} and \code{show} methods to make displaying large matrices better. \code{print} accepts an additional argument \code{n} corresponding to the maximum number of rows/columns to print before truncating.
\item \code{simMat} objects support matrix-style get/set operations like \code{s[1,]} or \code{s[1,3:5]}
\item \code{simMat} objects allow any values on the diagonal, rather than just
zeros as in \code{dist} objects.
\item \code{simMat} objects support conversion to matrices and \code{data.frame} objects
\item \code{simMat} objects implement get/set \code{Diag()} methods. Note usage of capitalized \code{Diag}; this is to avoid conflicts and weirdness with using base \code{diag}.
}

See the examples for details on using these features.
}
\value{
\code{simMat} and \code{as.simMat} return an object of class \code{"simMat"}. Internally,
the object stores the upper triangle of the matrix similar to how \code{dist}
stores objects. 

The object has the following attributes (besides \code{"class"} equal to \code{"simMat"}):

\item{nrow}{ the number of rows in the matrix implied by the vector }
\item{NAMES}{ the names of the rows/columns }

\code{as.matrix(s)} returns the equivalent matrix to a \code{"simMat"} object.

\code{as.data.frame(s)} returns a \code{data.frame} object corresponding to
pairwise similarities.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\examples{
## Creating a blank simMat object initialized to zeros
s <- simMat(0, nelem=20)
s

## Print out 5 rows instead of 10
print(s, n=5)

## Create a simMat object with 5 entries from a vector
dimn <- 5
vec <- 1:(dimn*(dimn-1) / 2)
s1 <- as.simMat(vec, DIAG=FALSE)
s1

## Here we include the diagonal
vec <- 1:(dimn*(dimn+1) / 2)
s2 <- as.simMat(vec, DIAG=TRUE)
s2

## Subsetting
s2[1,]
s2[1,3:4]
# all entries except first row
s2[-1,]
# all combos not including 1
s2[-1,-1]

## Replace values (automatically recycled)
s2[1,] <- 10
s2

## Get/set diagonal
Diag(s1)
Diag(s1) <- 5
s1

}
