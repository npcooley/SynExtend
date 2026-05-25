\name{EvaluatePairs}
\alias{EvaluatePairs}
\title{
Evaluating and Filtering Candidate Feature Pairs Against Decoy Alignments
}
\description{
A function for evaluating candidate genomic feature pairs in a \code{PairSummaries} object by comparing them against a set of decoy alignments. Decoys may be generated on the fly from a supplied database, or provided as a pre-built \code{PairSummaries} object. An optional evaluation method appends a model-derived or cluster-derived score to each pair. When \code{FDRCriteria} is supplied, pairs are filtered to a user-specified false discovery rate based on the ranking of a chosen column; otherwise the augmented object is returned unfiltered.
}
\usage{
EvaluatePairs(InputPairs,
              DataBase01,
              StaticDecoys,
              DecoyScalar = 0.5,
              EvaluationMethod = c("none",
                                   "kmeans",
                                   "glm",
                                   "lm"),
              Verbose = FALSE,
              FDRCriteria = c("Delta_Background" = 0.001),
              ...)
}
\arguments{
  \item{InputPairs}{
An object of class \code{PairSummaries}, typically produced by \code{SummarizePairs}. This is the set of candidate feature pairs to be evaluated.
}
  \item{DataBase01}{
Optional. Either a connection object to a \code{DECIPHER}-compatible SQLite database, or a character string giving the path to such a database on disk. When supplied and \code{StaticDecoys} is absent, decoy alignments are generated on the fly via \code{CreateDecoys} using attributes inherited from \code{InputPairs} (\code{GeneCalls}, \code{KmerSize}, \code{DefaultTranslationTable}). If a path string is supplied, the \code{RSQLite} package must be installed. At least one of \code{DataBase01} or \code{StaticDecoys} must be provided when \code{FDRCriteria} is not \code{NULL}.
}
  \item{StaticDecoys}{
Optional. A pre-built \code{PairSummaries} object to be used as the decoy set in place of on-the-fly generation. Takes precedence over \code{DataBase01} when both are supplied. Must be of class \code{PairSummaries}.
}
  \item{DecoyScalar}{
A numeric value strictly greater than zero and less than or equal to 1 controlling the size of the decoy set relative to \code{InputPairs}. The number of decoy rows sampled or generated is \code{ceiling(nrow(InputPairs) * DecoyScalar)}. When generating decoys on the fly, this ceiling is passed as \code{AA_limit} to \code{CreateDecoys}; if fewer decoys are available than requested, all available decoys are used. Defaults to \code{0.5}.
}
  \item{EvaluationMethod}{
A character string specifying the method used to append a \code{criteria_value} column to the combined candidate and decoy data. Must be one of:
\describe{
  \item{\code{"none"}}{No model is fitted and no \code{criteria_value} column is appended. The function returns the augmented \code{PairSummaries} object directly (after any FDR filtering). This is the default.}
  \item{\code{"glm"}}{A quasi-binomial generalised linear model is fitted with \code{Response} as the outcome and \code{Consensus}, \code{FeatureDiff}, \code{MaxMatch}, \code{KDist}, \code{Local_PID}, \code{Approx_Global_PID}, and \code{Delta_Background} as predictors. Fitted values from the model are stored in \code{criteria_value}.}
  \item{\code{"lm"}}{A linear model is fitted using the same predictors as \code{"glm"}, with \code{Delta_Background} as the response. Fitted values are stored in \code{criteria_value}.}
  \item{\code{"kmeans"}}{K-means clustering is applied across the predictor columns (without a \code{Response} column). The optimal number of clusters is selected by fitting a one-site binding curve to the within-cluster sum of squares elbow and scaling by \code{SelectScalar} (default 3, adjustable via \code{...}). Cluster assignments are stored in \code{criteria_value} and cluster centres are stored as the \code{centers} attribute of the returned object.}
}
Partial matching is supported.
}
  \item{Verbose}{
Logical indicating whether to print status messages and the elapsed time upon completion. Defaults to \code{FALSE}.
}
  \item{FDRCriteria}{
Either \code{NULL} or a named numeric vector of length 1 specifying the false discovery rate threshold to apply. The name must match a column of the combined evaluation data (e.g. \code{c("Delta_Background" = 0.001)}). Rows are ranked in decreasing order of the named column; the cumulative proportion of decoy rows in the ranking is computed, and all rows at or beyond the point where that proportion first reaches the threshold value are discarded. Only non-decoy rows in the retained set are returned. Requires decoy alignments to be present (either via \code{DataBase01} or \code{StaticDecoys}). Set to \code{NULL} to skip filtering. Defaults to \code{c("Delta_Background" = 0.001)}.
}
  \item{\dots}{
Additional named arguments passed to internal helpers. Currently recognised optional arguments include \code{K_val_02} (k-mer width for amino acid index construction in on-the-fly decoy generation; defaults to \code{6}) and \code{K_val_03} (k-mer width for nucleotide index construction; defaults to \code{10}). For \code{EvaluationMethod = "kmeans"}, \code{MaxK} (maximum number of clusters to test; default \code{15}) and \code{SelectScalar} (scaling factor applied to the estimated elbow; default \code{3}) are also accepted.
}
}
\details{
The function proceeds in three stages.

\strong{Stage 1 — Decoy assembly.} If \code{DataBase01} is supplied without \code{StaticDecoys}, decoy alignments are generated by calling \code{CreateDecoys} with arguments inherited from the attributes of \code{InputPairs}. If \code{StaticDecoys} is supplied, it is used directly. In both cases, the decoy set is downsampled to \code{ceiling(nrow(InputPairs) * DecoyScalar)} rows (or all available rows if fewer exist). Candidate pairs receive a \code{Response} value of \code{1}; decoy pairs receive \code{0}. If neither source of decoys is supplied, all rows are assigned \code{Response = 1} and the function proceeds without a null distribution.

\strong{Stage 2 — Evaluation.} A derived predictor \code{FeatureDiff} is computed as the absolute difference in feature lengths divided by the length of the longer feature, giving a normalised measure of length asymmetry. \code{MaxMatch} is similarly normalised by the sum of the two feature lengths. Depending on \code{EvaluationMethod}, a model or clustering procedure is applied to the combined data and a \code{criteria_value} column is appended.

\strong{Stage 3 — FDR filtering.} When \code{FDRCriteria} is not \code{NULL}, pairs are ranked in decreasing order of the named column and a cumulative FDR is computed as the running proportion of decoy rows. All rows at and beyond the rank position where the FDR first meets or exceeds the threshold are discarded. Decoy rows are removed from the result unconditionally; only candidate rows that survive the threshold are returned.
}
\value{
An object of class \code{c("data.frame", "PairSummaries")} containing the candidate pairs that survived evaluation and, if \code{FDRCriteria} was not \code{NULL}, FDR filtering. Columns include all columns from \code{InputPairs} plus:

\describe{
  \item{Response}{Integer. \code{1} for candidate pairs and \code{0} for decoy pairs. Decoy rows are always removed before the object is returned when \code{FDRCriteria} is not \code{NULL}; when \code{FDRCriteria} is \code{NULL} the combined candidate-and-decoy data is returned and this column is present.}
  \item{criteria_value}{Numeric or integer. Present only when \code{EvaluationMethod} is not \code{"none"}. Contains fitted values from the GLM or LM, or integer cluster assignments from k-means.}
}

The returned object retains the following attributes inherited from \code{InputPairs}:
\describe{
  \item{GeneCalls}{The named list of gene calls.}
  \item{KmerSize}{The k-mer size used in the upstream \code{SummarizePairs} call.}
  \item{DefaultTranslationTable}{The translation table identifier used upstream.}
}

When \code{EvaluationMethod = "kmeans"}, the returned object additionally carries a \code{centers} attribute containing the k-means cluster centre matrix.
}
\author{
  Nicholas Cooley \email{Nicholas.Cooley@ul.ie}
}
\seealso{
  \itemize{
    \item \code{\link{SummarizePairs}}
    \item \code{\link{CreateDecoys}}
    \item \code{\link{NucleotideOverlap}}
  }
}
\examples{
library(DBI)
data("init_pairs")
tmp01 <- system.file("extdata",
                      "example_db.sqlite",
                      package = "SynExtend")
tmp02 <- tempfile()
file.copy(from = tmp01,
          to = tmp02)

drv <- dbDriver("SQLite")
conn01 <- dbConnect(drv = drv,
                    tmp02)
                    
x <- EvaluatePairs(InputPairs = init_pairs,
                   DataBase01 = conn01)
}
