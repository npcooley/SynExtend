\name{predict.EvoWeaver}
\alias{predict.EvoWeaver}
\title{
Make predictions with EvoWeaver objects
}
\description{
This S3 method predicts pairwise functional associations between gene groups encoded in a \code{\link{EvoWeaver}} object.
This returns an object of type \code{\link{EvoWeb}}, which is essentially an adjacency
matrix with some extra S3 methods to make printing cleaner.
}
\usage{
\method{predict}{EvoWeaver}(object, Method='Ensemble',
         Subset=NULL, Processors=1L,
         MySpeciesTree=SpeciesTree(object, Verbose=Verbose),
         PretrainedModel="KEGG",
         NoPrediction=FALSE,
         ReturnDataFrame=TRUE,
         Verbose=interactive(), ...)
}
\arguments{
  \item{object}{
    A EvoWeaver object
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

    \code{subset=rbind(c(1,2),c(1,3),c(2,3))} produces equivalent functionality to
    \code{subset=1:3}.
  }
  \item{Processors}{
    Number of cores to use for methods that support multithreaded execution.
    Setting to \code{NULL} or a negative value will use the value of
    \code{detectCores()},
    or one core if the number of available cores cannot be determined. See Note
    for more information.
  }
  \item{MySpeciesTree}{
    Phylogenetic tree of all genomes in the dataset. Required for \code{Method=c('RPContextTree', 'GLDistance', 'CorrGL', 'MoransI', 'Behdenna')}. \code{'Behdenna'} requires a rooted, bifurcating tree (other values of \code{Method} can handle arbitrary trees). Note that \code{EvoWeaver} can automatically infer a species tree if initialized with \code{dendrogram} objects.
  }
  \item{PretrainedModel}{
    A pretrained model for use with ensemble predictions. The default
    value is \code{"KEGG"}, corresponding to a built-in ensemble model
    trained on the KEGG MODULE database. Alternative values allowed
    are \code{"CORUM"}, for a built-in ensemble model trained on the
    CORUM database, or any user-trained model.
    See the examples for how to train an
    ensemble method to pass to \code{PretrainedModel}.

    Has no effect if \code{Method != 'Ensemble'}.
  }
  \item{NoPrediction}{
    For \code{Method='Ensemble'}, should data be returned prior to making predictions?

    If \code{TRUE}, this will instead return a \link[base]{data.frame} object
    with predictions from each algorithm for each pair. This dataframe is typically
    used to train an ensemble model.

    If \code{FALSE}, EvoWeaver will return predictions for each pair (using user model
    if provided or a built-in otherwise).
  }
  \item{ReturnDataFrame}{
    Logical indicating whether to return a \code{data.frame} object or a list of \code{EvoWeb} objects.
    Defaults to \code{TRUE}. Setting this parameter to \code{FALSE} is not recommended for typical users.
  }
  \item{Verbose}{
    Logical indicating whether to print progress bars and messages. Defaults to \code{TRUE}.
  }
  \item{...}{
    Additional parameters for other predictors and consistency with generic.
  }
}

\details{
\code{predict.EvoWeaver} wraps several methods to create an easy interface for multiple prediction types. \code{Method='Ensemble'} is the default value, but each of the component analyses can also be accessed. Common arguments to \code{Method} include:

\itemize{
  \item \code{'Ensemble'}: Ensemble prediction combining individual coevolutionary predictors. See \code{Note} below.
  \item \code{'PhylogeneticProfiling'}: All \link[=EvoWeaver-PPPreds]{Phylogenetic Profiling Algorithms} used in the EvoWeaver manuscript.
  \item \code{'PhylogeneticStructure'}: All \link[=EvoWeaver-PSPreds]{EvoWeaver Phylogenetic Structure Methods}
  \item \code{'GeneOrganization'}: All \link[=EvoWeaver-PSPreds]{EvoWeaver Gene Organization Methods}
  \item \code{'SequenceLevel'}: All \link[=EvoWeaver-SLPreds]{EvoWeaver Sequence Level Methods} used in the EvoWeaver manuscript.
}

Additional information and references for each prediction algorithm can be found at the
following pages:

\itemize{
\item \link[=EvoWeaver-PPPreds]{EvoWeaver Phylogenetic Profiling Methods}

\item \link[=EvoWeaver-PSPreds]{EvoWeaver Phylogenetic Structure Methods}

\item \link[=EvoWeaver-GOPreds]{EvoWeaver Gene Organization Methods}

\item \link[=EvoWeaver-SLPreds]{EvoWeaver Sequence-Level Methods}
}

The standard return type is a \code{data.frame} object with one column per predictor and an additional two columns specifying the genes in each pair. If \code{ReturnDataFrame=FALSE}, this returns a \code{EvoWeb} object. See \code{\link{EvoWeb}} for more information. Use of this parameter is discouraged.

Different methods require different types of input. The constructor
\code{\link{EvoWeaver}} will notify the user which methods are
runnable with the given data. Method \code{Ensemble} automatically selects the
methods that can be run with the given input data.

See \code{\link{EvoWeaver}} for more information on input data types.

Complete listing of all supported methods (asterisk denotes a method used in \code{Ensemble}, if possible):
\itemize{
  \item \code{'ExtantJaccard'}: Jaccard Index of Presence/Absence (P/A) profiles at extant leaves
  \item \code{'Hamming'}: Hamming similarity of P/A profiles
  \item * \code{'GLMI'}: MI of G/L profiles
  \item \code{'PAPV'}: \code{1-p_value} of P/A profiles
  \item \code{'ProfDCA'}: Direct Coupling Analysis of P/A profiles
  \item \code{'Behdenna'}: Analysis of Gain/Loss events following Behdenna et al. (2016)
  \item \code{'CorrGL'}: Correlation of ancestral Gain/Loss events
  \item * \code{'GLDistance'}: Score-based method based on distance between inferred ancestral Gain/Loss events
  \item * \code{'PAJaccard'}: Centered Jaccard distance of P/A profiles with conserved clades collapsed
  \item * \code{'PAOverlap'}: Conservation of ancestral states based on P/A profiles
  \item * \code{'RPMirrorTree'}: MirrorTree using Random Projection for dimensionality reduction
  \item * \code{'RPContextTree'}: MirrorTree with Random Projection correcting for species tree and P/A conservation
  \item * \code{'GeneDistance'}: Co-localization analysis
  \item * \code{'MoransI'}: Co-localization analysis using   \link[=MoranI]{Moran's I} for phylogenetic correction and significance
  \item * \code{'OrientationMI'}: Mutual Information of Gene Relative Orientation
  \item * \code{'GeneVector'}: Correlation of distribution of sequence level residues following Zhao et al. (2022)
  \item * \code{'SequenceInfo'}: Mutual information of sites in multiple sequence alignment
}
}
\value{
Returns a \code{data.frame} object where each row corresponds to a single prediction for a pair of gene groups. The first two columns contain the gene group identifiers for each pair, and the remaining columns contain each prediction.

If \code{ReturnDataFrame=FALSE}, the return type is a list of \code{EvoWeb} objects. See \code{\link{EvoWeb}} for more info.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\note{
EvoWeaver's publication used a random forest model from the \code{randomForest} package for prediction. The next release of EvoWeaver will include multiple new built-in ensemble methods, but in the interim users are recommended to rely on \code{randomForest} or \code{neuralnet}. Planned algorithms are random forests and feed-forward neural networks. Feel free to contact me regarding other models you would like to see added.

  If \code{NumCores} is set to \code{NULL}, EvoWeaver will use one less core than is detected, or one core if \code{detectCores()}
  cannot detect the number of available cores. This is because of a recurring issue
  on my machine where the R session takes all available cores and is then locked
  out of forking processes, with the only solution to restart the entire R session.
  This may be an issue specific to ARM Macs, but out of an abundance of caution
  I've made the default setting to be slightly slower but guarantee completion
  rather than risk bricking a machine.
}

\seealso{
\code{\link{EvoWeaver}}

\code{\link{EvoWeb}}

\link[=EvoWeaver-PPPreds]{EvoWeaver Phylogenetic Profiling Predictors}

\link[=EvoWeaver-PSPreds]{EvoWeaver Phylogenetic Structure Predictors}

\link[=EvoWeaver-GOPreds]{EvoWeaver Gene Organization Predictors}

\link[=EvoWeaver-SLPreds]{EvoWeaver Sequence-Level Predictors}
}
\examples{
##############
## Prediction with built-in model and data
###############

set.seed(555L)
exData <- get(data("ExampleStreptomycesData"))
ew <- EvoWeaver(exData$Genes[1:50], MySpeciesTree=exData$Tree)

# Subset isn't necessary but is faster for a working example
evoweb1 <- predict(ew, Subset=1:2)

# print out results as an adjacency matrix
if(interactive()) print(evoweb1)

###############
## Training own ensemble model
###############

datavals <- evoweb1[,-c(1,2,10)]
actual_values <- sample(c(0,1), nrow(datavals), replace=TRUE)
# This example just picks random numbers
# ***Do not do this for your own models***

# Make sure the actual values correspond to the right pairs!
datavals[,'y'] <- actual_values
myModel <- glm(y~., datavals[,-c(1,2)], family='binomial')

testEvoWeaverObject <- EvoWeaver(exData$Genes[51:60], MySpeciesTree=exData$Tree)
evoweb2 <- predict(testEvoWeaverObject,
                     PretrainedModel=myModel)

# Print result as a data.frame of pairwise scores
if(interactive()) print(evoweb2)
}
