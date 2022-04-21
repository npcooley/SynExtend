\name{plot.ProtWeb}
\alias{plot.ProtWeb}
\title{
Plot predictions in a ProtWeb object
}
\description{
ProtWeb objects are outputted from \code{\link{predict.ProtWeaver}}. 

This function plots the predictions in the object using a force-directed embedding
of connections in the adjacency matrix.

\emph{This function is still a work in progress.}
}
\usage{
\method{plot}{ProtWeb}(x, nsim=10, gravity=0.05, 
      coulomb=0.1, connection=5,
      move_rate=0.25, cutoff=0.2, 
      verbose=TRUE, colorpalette=topo.colors, ...)
}
\arguments{
  \item{x}{
    A ProtWeb object. See \code{\link{predict.ProtWeaver}}
  }
  \item{nsim}{
    Number of iterations to run the model for.
  }
  \item{gravity}{
    Strength of Gravity force. See 'Details'.
  }
  \item{coulomb}{
    Strength of Coulomb force. See 'Details'.
  }
  \item{connection}{
    Strength of Connective force. See 'Details'.
  }
  \item{move_rate}{
    Controls how far each point moves in each iteration.
  }
  \item{cutoff}{
    Cutoff value; if \code{abs(val) < cutoff}, that connection is shrunk to zero.
  }
  \item{verbose}{
    Should progress be printed?
  }
  \item{colorpalette}{
    Color palette for graphing. Valid inputs are any palette available in
    \code{palette.pals()}. See \code{\link[grDevices]{palette}} for more info.
  }
  \item{...}{
    Additional parameters for consistency with generic.
  }
}

\details{
  This function plots the \code{ProtWeb} object using a force-directed embedding.
  This embedding has three force components:
  \itemize{
  \item Gravity Force: Attractive force pulling nodes towards \code{(0,0)}
  \item Coulomb Force: Repulsive force pushing close nodes away from each other
  \item Connective Force: Tries to push node connections to equal corresponding values in the adjacency matrix
  }
  
  The parameters in the function are sufficient to get an embedding, though users
  are welcome to try to tune them for a better visualization. This function is meant
  to aid with visualization of the adjacency matrix, not for concrete analyses of clusters.
  
  I'm not satisfied with the speed or visualization of this embedding, so next release
  cycle will update this function with an updated version of this algorithm to improve
  plotting.
}
\value{
No return value; creates a plot in the graphics window.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\seealso{
\code{\link{predict.ProtWeaver}}

\code{\link{GetProtWebData.ProtWeb}}
}
\examples{
exData <- get(data("ExampleStreptomycesData"))
pw <- ProtWeaver(exData$Genes)

# Subset isn't necessary but is faster for a working example
# Same w/ method='Jaccard'
protweb <- predict(pw, 'Jaccard', subset=1:50)

# plot(protweb)
}