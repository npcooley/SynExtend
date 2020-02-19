\name{ExactSelect}
\alias{ExactSelect}
\docType{data}
\title{
Model for identifying erroneously linked pairs
}
\description{
Though the function PairSummaries provides an argument allowing users to ask for alignments, given the time consuming nature of that process on large datasets, models are provided which allow for the quick and efficient identification of pairs whose PID would likely fall within a random distribution of PIDs.
}
\usage{data("ExactSelect")}
\format{
  The format is:
List of 21
 $ coefficients : Named num [1:32] 0.662363 46.993855 -2.187871 -0.000554 -4.981938 ...
  ..- attr(*, "names")= chr [1:32] "(Intercept)" "MaxCoverage" "MinCoverage" "ExactMatch" ...
 $ R            : num [1:32, 1:32] -63.5 0 0 0 0 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:32] "(Intercept)" "MaxCoverage" "MinCoverage" "ExactMatch" ...
  .. ..$ : chr [1:32] "(Intercept)" "MaxCoverage" "MinCoverage" "ExactMatch" ...
 $ rank         : int 32
 $ qr           :List of 4
  ..$ rank : int 32
  ..$ qraux: num [1:32] 1 1 1 1 1 ...
  ..$ pivot: int [1:32] 1 2 3 4 5 6 7 8 9 10 ...
  ..$ tol  : num 1e-11
  ..- attr(*, "class")= chr "qr"
 $ family       :List of 7
  ..$ family    : chr "binomial"
  ..$ link      : chr "logit"
  ..$ linkfun   :function (mu)  
  ..$ linkinv   :function (eta)  
  ..$ mu.eta    :function (eta)  
  ..$ initialize:  expression({  if (NCOL(y) == 1) {  if (is.factor(y))  y <- y != levels(y)[1L]  n <- rep.int(1, nobs)  y[weights =| __truncated__}}
  ..$ valideta  :function (eta)  
  ..- attr(*, "class")= chr "family"
 $ deviance     : num 31409
 $ aic          : num 31473
 $ null.deviance: num 72189
 $ iter         : int 9
 $ df.residual  : int 196744
 $ df.null      : int 196775
 $ converged    : logi TRUE
 $ boundary     : logi FALSE
 $ call         : language glm(formula = as.formula(ExactModels[[1]][[31]]), family = binomial(link = "logit"), data = TrainingTable02)
 $ formula      :Class 'formula'  language ExactCode ~ MaxCoverage * MinCoverage * ExactMatch * NormDeltaStart * NormDeltaStop
 $ terms        :Classes 'terms', 'formula'  language ExactCode ~ MaxCoverage * MinCoverage * ExactMatch * NormDeltaStart * NormDeltaStop
  .. ..- attr(*, "variables")= language list(ExactCode, MaxCoverage, MinCoverage, ExactMatch, NormDeltaStart, NormDeltaStop)
  .. ..- attr(*, "factors")= int [1:6, 1:31] 0 1 0 0 0 0 0 0 1 0 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:6] "ExactCode" "MaxCoverage" "MinCoverage" "ExactMatch" ...
  .. .. .. ..$ : chr [1:31] "MaxCoverage" "MinCoverage" "ExactMatch" "NormDeltaStart" ...
  .. ..- attr(*, "term.labels")= chr [1:31] "MaxCoverage" "MinCoverage" "ExactMatch" "NormDeltaStart" ...
  .. ..- attr(*, "order")= int [1:31] 1 1 1 1 1 2 2 2 2 2 ...
  .. ..- attr(*, "intercept")= int 1
  .. ..- attr(*, "response")= int 1
  .. ..- attr(*, "predvars")= language list(ExactCode, MaxCoverage, MinCoverage, ExactMatch, NormDeltaStart, NormDeltaStop)
  .. ..- attr(*, "dataClasses")= Named chr [1:6] "numeric" "numeric" "numeric" "numeric" ...
  .. .. ..- attr(*, "names")= chr [1:6] "ExactCode" "MaxCoverage" "MinCoverage" "ExactMatch" ...
 $ offset       : NULL
 $ control      :List of 3
  ..$ epsilon: num 1e-08
  ..$ maxit  : num 25
  ..$ trace  : logi FALSE
 $ method       : chr "glm.fit"
 $ contrasts    : NULL
 $ xlevels      : Named list()
 - attr(*, "class")= chr [1:2] "glm" "lm"
}
\details{
A model for rejecting identified pairs whose link statistics indicate a likely exact PID that would fall within a random distribution in AA Space.
}
\source{
TBA
}
\references{
TBA
}
\examples{
data(ExactSelect)
## ? ?
}
\keyword{datasets}
