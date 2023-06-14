\name{predict.ProtWeaver}
\alias{predict.ProtWeaver}
\title{
Make predictions with ProtWeaver objects
}
\description{
This S3 method predicts pairwise functional associations between gene groups encoded in a \code{\link{ProtWeaver}} object.
This returns an object of type \code{\link{ProtWeb}}, which is essentially an adjacency
matrix with some extra S3 methods to make printing cleaner.
}
\usage{
\method{predict}{ProtWeaver}(object, Method='Ensemble',
         Subset=NULL, Processors=1L,
         MySpeciesTree=SpeciesTree(object), 
         PretrainedModel=NULL,
         NoPrediction=FALSE,
         ReturnRawData=FALSE, Verbose=TRUE, ...)
}
\arguments{
  \item{object}{
    A ProtWeaver object
  }
  \item{Method}{
    Method(s) to use for prediction. This can be a character vector with multiple entries for predicting using multiple methods. See 'Details' for more information.
  }
  \item{Subset}{
    Subset of data to predict on. This can either be a vector or a \code{2xN} matrix.
    
    
    If a vector, prediction proceeds for all possible pairs of elements specified in the vector 
    (either by name, for character vector, or by index, for numeric vector). For example, 
    \code{subset=1:3} will predict for pairs \code{(1,2), (1,3), (2,3)}. 
    
    If a matrix, subset is interpreted as a matrix of pairs, where each row of the matrix 
    specifies a pair to evaluate. These can also be specifed by name (character) or by index 
    (numeric). 
    
    \code{subset=cbind(c(1,1,2), c(2,3,3))} produces equivalent functionality to 
    \code{subset=1:3}.
  }
  \item{Processors}{
    Number of cores to use for methods that support multithreaded execution. 
    Currently only supported for methods \code{'ProfDCA'}, \code{'MirrorTree'} and \code{'Ensemble'}. Setting
    to \code{NULL} or a negative value will use the value of \code{detectCores()}, 
    or one core if the number of available cores cannot be determined. See Note
    for more information. This parameter has no effect on Windows due to reliance on
    forking via \code{mclapply}.
  }
  \item{MySpeciesTree}{
    Phylogenetic tree of all genomes in the dataset. Required for \code{Method='Behdenna'},
    and can improve predictions for other methods. \code{'Behdenna'} requires a 
    rooted, bifurcating tree (other values of \code{Method} can handle arbitrary trees). Note that \code{ProtWeaver} can automatically infer a species tree if initialized with \code{dendrogram} objects.
  }
  \item{PretrainedModel}{
    A pretrained model for use with ensemble predictions. If unspecified when 
    \code{Method='Ensemble'}, the program will use built-in models 
    (see \code{\link{BuiltInEnsembles}}). See the examples for how to train an
    ensemble method to pass to \code{PretrainedModel}. 
    
    Has no effect if \code{Method != 'Ensemble'}.
  }
  \item{NoPrediction}{
    For \code{Method='Ensemble'}, should data be returned prior to making predictions?
    
    If \code{TRUE}, this will instead return a \link[base]{data.frame} object
    with predictions from each algorithm for each pair. This dataframe is typically
    used to train an ensemble model. 
    
    If \code{FALSE}, ProtWeaver will return predictions for each pair (using user model
    if provided or a built-in otherwise).
  }
  \item{ReturnRawData}{
      Internal parameter used for ensemble predictions. Should not be set by the user.
  }
  \item{Verbose}{
    Logical indicating whether to print progress bars and messages. Defaults to \code{TRUE}.
  }
  \item{...}{
    Additional parameters for other predictors and consistency with generic.
  }
}

\details{
\code{predict.ProtWeaver} wraps several methods to create an easy interface for multiple prediction types. \code{Method='Ensemble'} is the default value, but the following values of \code{Method} are also supported: 

\itemize{
  \item \code{'Ensemble'}: Ensemble prediction combining individual coevolutionary predictors.
  \item \code{'Jaccard'}: Jaccard distance of Presence/Absence (P/A) profiles
  \item \code{'Hamming'}: Hamming distance of P/A profiles
  \item \code{'MutualInformation'}: MI of P/A profiles
  \item \code{'PAPV'}: \code{1-p_value} of P/A profiles
  \item \code{'ProfDCA'}: Direct Coupling Analysis of P/A profiles
  \item \code{'Behdenna'}: Analysis of Gain/Loss events following Behdenna et al. (2016) 
  \item \code{'CorrGL'}: Correlation of ancestral Gain/Loss events
  \item \code{'GainLoss'}: Score-based method based on distance between inferred ancestral Gain/Loss events
  \item \code{'MirrorTree'}: MirrorTree using Random Projection for dimensionality reduction
  \item \code{'ContextTree'}: MirrorTree with Random Projection correcting for species tree and P/A conservation
  \item \code{'Coloc'}: Co-localization analysis
  \item \code{'ColocMoran'}: Co-localization analysis using   \link[=MoransI]{Moran's I} for phylogenetic correction and significance
  \item \code{'TranscripMI'}: Mutual Information of Transcriptional Direction
  \item \code{'NVDT'}: Correlation of distribution of sequence level residues following Zhao et al. (2022)
}

The best performing individual predictors are \code{c('CorrGL', 'GainLoss', 'MirrorTree', 'Jaccard')}. Users interesting in running quick analyses should use \code{c('CorrGL', 'GainLoss', 'Jaccard')}.  

Additional information and references for each prediction algorithm can be found at the 
following pages:

\itemize{
\item \link[=ProtWeaver-PAPreds]{ProtWeaver Presence/Absence Methods}

\item \link[=ProtWeaver-DMPreds]{ProtWeaver Distance Matrix Methods}

\item \link[=ProtWeaver-ColocPreds]{ProtWeaver Co-localization Methods}

\item \link[=ProtWeaver-ResiduePreds]{ProtWeaver Residue Level Methods}
}

This returns a \code{ProtWeb} object, an S3 class that makes formatting and printing of 
results slightly nicer. See \code{\link{ProtWeb}} for more information.

Different methods require different types of input. The constructor 
\code{\link{ProtWeaver}} will notify the user which methods are
runnable with the given data. Method \code{Ensemble} automatically selects the 
methods that can be run with the given input data. 

See \code{\link{ProtWeaver}} for more information on input data types.

}
\value{
Returns a ProtWeb object. See \code{\link{ProtWeb}} for more info.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\note{
  \code{NumCores} uses 1 less core than is detected, or 1 core if \code{detectCores()}
  cannot detect the number of available cores. This is because of a recurring issue
  on my machine where the R session takes all available cores and is then locked
  out of forking processes, with the only solution to restart the entire R session.
  This may be an issue specific to ARM Macs, but out of an abundance of caution
  I've made the default setting to be slightly slower but guarantee completion
  rather than risk bricking a machine.

  More models will be implemented in the future. Planned models for next release
  include:
  \itemize{
    \item Random Forests for Ensemble predictions
    \item XGBoost for Ensemble predictions
  }
  
  Feel free to contact me regarding other models you would like to see added.
}

\seealso{
\code{\link{ProtWeaver}}

\code{\link{ProtWeb}}

\link[=ProtWeaver-PAPreds]{ProtWeaver Presence/Absence Predictors}

\link[=ProtWeaver-DMPreds]{ProtWeaver Distance Matrix Predictors}

\link[=ProtWeaver-ColocPreds]{ProtWeaver Co-localization Predictors}

\link[=ProtWeaver-ResiduePreds]{ProtWeaver Residue Level Predictors}
}
\examples{
##############
## Prediction with built-in model and data 
###############

exData <- get(data("ExampleStreptomycesData"))
pw <- ProtWeaver(exData$Genes[1:50])

# Subset isn't necessary but is faster for a working example
protweb1 <- predict(pw, Subset=1:10, MySpeciesTree=exData$Tree)

# print out results as an adjacency matrix
protweb1

###############
## Training own ensemble model
###############

datavals <- predict(pw, NoPrediction=TRUE)
                  
actual_values <- sample(c(0,1), nrow(datavals), replace=TRUE)
# This example just picks random numbers
# ***Do not do this for your own models***

# Make sure the actual values correspond to the right pairs! 
datavals[,'y'] <- actual_values                  
myModel <- glm(y~., datavals[,-c(1,2)], family='binomial')

testProtWeaverObject <- ProtWeaver(exData$Genes[51:60])
protweb2 <- predict(testProtWeaverObject, 
                     PretrainedModel=myModel)
                     
# Print result as a matrix of pairwise scores
protweb2
}