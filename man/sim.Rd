\name{sim}
\alias{sim}
\alias{sim-class}
\alias{as.sim}
\alias{diag.sim}
\title{
Similarity Matrices
}
\description{
The \code{sim} object is an internally utilized class that provides similar 
functionality to the \code{\link{dist}} object, but with matrix-like accessors.

Like \code{dist}, this object stores values as a vector, reducing memory by
making use of assumed symmetry. \code{sim} currently only supports numeric
data types.
}
\usage{
## Create a blank sym object
sym(VALUE, nelem, NAMES=NULL, DIAG=FALSE)

## S3 methods for 'vector' and 'matrix':
as.sim(m, NAMES=NULL, DIAG=TRUE)

## S3 method for class 'sim'
print(s, n=10)

## S3 method for class 'sim'
as.matrix(s)
as.data.frame(s)

## S3 method for class 'sim'
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
}
\details{
The \code{sim} object has a very similar format to \code{dist} objects, but 
with a few notable changes:

\itemize{
\item \code{sim} objects have streamlined \code{print} and \code{show} methods to make displaying large matrices better. \code{print} accepts an additional argument \code{n} corresponding to the maximum number of rows/columns to print before truncating.
\item \code{sim} objects support matrix-style get/set operations like \code{s[1,]} or \code{s[1,3:5]}
\item \code{sim} objects allow any values on the diagonal, rather than just
zeros as in \code{dist} objects.
\item \code{sim} objects support conversion to matrices and \code{data.frame} objects
\item \code{sim} objects implement get/set \code{diag()} methods 
}

See the examples for details on using these features.
}
\value{
\code{sim} and \code{as.sim} return an object of class \code{"sim"}. Internally,
the object stores the upper triangle of the matrix similar to how \code{dist}
stores objects. 

The object has the following attributes (besides \code{"class"} equal to \code{"sim"}):

\item{nrow}{ the number of rows in the matrix implied by the vector }
\item{NAMES}{ the names of the rows/columns }

\code{as.matrix(s)} returns the equivalent matrix to a \code{"sim"} object.

\code{as.data.frame(s)} returns a \code{data.frame} object corresponding to
pairwise similarities.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\examples{
## Creating a blank sim object initialized to zeros
s <- sim(0, nelem=20)
s

## Print out 5 rows instead of 10
print(s, n=5)

## Create a sim object with 5 entries from a vector
dimn <- 5
vec <- 1:(dimn*(dimn-1) / 2)
s1 <- as.sim(vec, DIAG=FALSE)
s1

## Here we include the diagonal
vec <- 1:(dimn*(dimn+1) / 2)
s2 <- as.sim(vec, DIAG=TRUE)
s2

## Subsetting
s2[1,]
s2[1,3:4]

## Replace values (automatically recycled)
s2[1,] <- 10
s2

## Get diagonal
diag(s1)
diag(s1) <- 5
s1

}
