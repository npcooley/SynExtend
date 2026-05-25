###### -- Given a pair summaries object, evaluate the pairs present -----------

###### -- NOTES ---------------------------------------------------------------
# this function technically only evaluates when 'FDRCriteria' is not NULL
# when it is not, and decoy alignments are supplied through StaticDecoys, or
# generated from a supplied 'DataBase01', if FDRCriteria is NULL it simply returns
# the input pairs with an appended column 'criteria_val' associated with the method

EvaluatePairs <- function(InputPairs,
                          DataBase01, # create decoys on the fly
                          StaticDecoys, # use a supplied SummarizePairs object as decoys
                          DecoyScalar = 0.5, # only sets limits based on AA rows
                          EvaluationMethod = c("none",
                                               "kmeans",
                                               "glm",
                                               "lm"),
                          Verbose = FALSE,
                          FDRCriteria = c("Delta_Background" = 0.001),
                          ...) {
  
  # requires decoy pairs
  append_by_glm <- function(Candidates,
                            SelectColumns = c("Consensus",
                                              "FeatureDiff",
                                              "MaxMatch",
                                              "KDist",
                                              "Local_PID",
                                              "Approx_Global_PID",
                                              "Delta_Background",
                                              "Response")) {
    # predict response based on the other columns
    # requires more than one response value
    Candidates <- Candidates[, SelectColumns]
    resp <- glm(formula = Response ~ .,
                family = "quasibinomial",
                data = Candidates)
    return(resp)
    # return a single column data frame of predicted response values
  }
  # predict the background score
  append_by_lm <- function(Candidates,
                           SelectColumns = c("Consensus",
                                             "FeatureDiff",
                                             "MaxMatch",
                                             "KDist",
                                             "Local_PID",
                                             "Approx_Global_PID",
                                             "Delta_Background"),
                           TargetColumn = "Delta_Background") {
    # does not require decoys
    Candidates <- Candidates[, SelectColumns]
    if (!(TargetColumn %in% SelectColumns)) {
      stop("Modeled value must be a column of the input data.")
    }
    formula_string <- paste(TargetColumn,
                            "~ .")
    resp <- lm(formula = as.formula(formula_string),
               data = Candidates)
    return(resp)
    # return a single column data.frame of the sorted predicted value
  }
  # does not expect a response column
  append_by_k <- function(Candidates,
                          SelectColumns = c("Consensus",
                                            "FeatureDiff",
                                            "MaxMatch",
                                            "KDist",
                                            "Local_PID",
                                            "Approx_Global_PID",
                                            "Delta_Background"),
                          MaxK = 15,
                          SelectScalar = 3,
                          NormCols = FALSE) {
    Candidates <- Candidates[, SelectColumns]
    
    if (NormCols) {
      # names should be retained here?
      Candidates <- as.data.frame(lapply(X = Candidates,
                                         FUN = function(x) {
                                           NormVec(x)
                                         }))
    }
    
    nclust <- seq(from = 2,
                  by = 1,
                  to = MaxK)
    kmc <- vector(mode = "list",
                  length = length(nclust))
    
    for (d1 in seq_along(kmc)) {
      kmc[[d1]] <- suppressWarnings(kmeans(x = Candidates,
                                           centers = nclust[d1],
                                           iter.max = 25L,
                                           nstart = 25L))
    }
    wss <- vapply(X = kmc,
                  FUN = function(x) {
                    x$tot.withinss
                  },
                  FUN.VALUE = vector(mode = "numeric",
                                     length = 1L))
    dat1 <- cbind("n" = nclust,
                  "wss" = wss)
    dat2 <- cbind("n" = nclust - 2,
                  "wss" = abs(wss - wss[1]))
    fita <- nls(dat2[, 2L]~OneSite(X = dat2[, 1L],
                                   Bmax,
                                   Kd),
                start = list(Bmax = max(dat2[, 2L]),
                             Kd = unname(quantile(dat2[, 1L], .25))))
    fitasum <- summary(fita)
    
    found_k <- ceiling((fitasum$coefficients["Kd", "Estimate"] + 1L) * SelectScalar)
    
    if (found_k > MaxK) {
      warning("Selected cluster number appears to be larger than the max number of clusters tested. Returning to max searched clusters.")
      found_k <- MaxK
    }
    if (found_k < 1L) {
      warning("Scalar selection requested a number of clusters less than 2, defaulting to 2 clusters.")
      found_k <- 1L
    }
    
    kmcselect <- kmc[[found_k]]
    # from here, i need to return this piece of data, and the rest of the 
    # gymnastics occur outside of this function...
    return(kmcselect)
  }
  
  # overhead checks for args
  if (Verbose) {
    tstart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  
  fun_dots <- list(...)
  
  if (DecoyScalar > 1 | DecoyScalar <= 0) {
    stop ("'DecoyScalar' must be less than 1 and greater than zero.")
  }
  
  if (!missingArg(StaticDecoys)) {
    if (!is(object = StaticDecoys,
            class2 = "PairSummaries")) {
      stop("Objects supplied to 'StaticDecoys' must be of class 'PairSummaries'.")
    }
  }
  
  # input pairs
  if (!is(object = InputPairs,
          class2 = "PairSummaries")) {
    stop ("'InputPairs' must be a object of class 'PairSummaries'.")
  }
  
  if (!is.null(FDRCriteria)) {
    if (length(FDRCriteria) > 1) {
      stop("If 'FDRCriteria' is not NULL it must be a named numeric vector of length 1.")
    }
  }
  
  # assign some things
  # these aren't needed because create decoys is doing this work now...
  # AA_matrix <- DECIPHER:::.getSubMatrix("PFASUM50")
  # NT_matrix <- DECIPHER:::.nucleotideSubstitutionMatrix(2L, -1L, 1L)
  
  # eval methods
  EvaluationMethod <- match.arg(EvaluationMethod)
  
  if (!missingArg(DataBase01) &
      missingArg(StaticDecoys)) {
    # if a database is supplied
    # and 
    # a static decoy set is not supplied
    # just build a generic decoy set
    if (is.character(DataBase01)) {
      if (!requireNamespace(package = "RSQLite",
                            quietly = TRUE)) {
        stop("Package 'RSQLite' must be installed.")
      }
      if (!("package:RSQLite" %in% search())) {
        print("Eventually character vector access to DECIPHER DBs will be deprecated.")
        requireNamespace(package = "RSQLite",
                         quietly = TRUE)
      }
      dbConn <- dbConnect(dbDriver("SQLite"), DataBase01)
      on.exit(dbDisconnect(dbConn))
      
    } else {
      dbConn <- DataBase01
      if (!dbIsValid(dbConn)) {
        stop("The connection has expired.")
      }
    }
    
    if (Verbose) {
      cat("\nGenerating decoy alignments.\n")
    }
    # CreateDecoys has some args set to NULL because I do unusual things,where
    # we can inherit some of the vals for those arguments from the attributes
    # of the input pairs object
    # we need to set the others
    # if (!("K_val_01" %in% fun_dots)) {
    #   fun_dots[["K_val_01"]] <- attr(x = InputPairs,
    #                                  which = "KmerSize")
    # }
    # there's no equivalent k size args to inherit from the input pairs
    # because we're not using searchindex in the initial search scheme
    if (!("K_val_02" %in% fun_dots)) {
      fun_dots[["K_val_02"]] <- 6
    }
    if (!("K_val_03" %in% fun_dots)) {
      fun_dots[["K_val_03"]] <- 10
    }
    fun_dots[["AA_limit"]] <- nrow(InputPairs) * DecoyScalar
    # create decoys doesn't subset, it just ends the search based on the limit,
    # so we still need to subset afterwards
    InputDecoys <- CreateDecoys(DataBase01 = dbConn,
                                GeneCalls = attr(x = InputPairs,
                                                 which = "GeneCalls"),
                                DefaultTranslationTable = attr(x = InputPairs,
                                                               which = "DefaultTranslationTable"),
                                K_val_01 = attr(x = InputPairs,
                                                which = "KmerSize"),
                                K_val_02 = fun_dots$K_val_02,
                                K_val_03 = fun_dots$K_val_03,
                                AA_limit = fun_dots$AA_limit,
                                Verbose = FALSE)
    curr_sample <- ceiling(nrow(InputPairs) * DecoyScalar)
    if (curr_sample > nrow(InputDecoys)) {
      curr_sample <- seq(nrow(InputDecoys))
    } else {
      curr_sample <- sample(x = nrow(InputDecoys),
                            size = curr_sample,
                            replace = FALSE)
    }
    InputDecoys <- InputDecoys[curr_sample, ]
    rownames(InputDecoys) <- NULL
    InputDecoys <- cbind(InputDecoys,
                         "Response" = rep(0, nrow(InputDecoys)))
    InputPairs <- cbind(InputPairs,
                        "Response" = rep(1, nrow(InputPairs)))
    EvalData <- rbind(InputPairs,
                      InputDecoys)
    
    
  } else if (!missingArg(StaticDecoys)) {
    # use a supplied set of static decoys as 'TRUE NEGATIVES'
    # must be a named list of length two with one position named 'AA' and the
    # other named 'NT'
    if (!is(object = StaticDecoys,
            class2 = "PairSummaries")) {
      stop("'StaticDecoys' must be an object of class 'PairSummaries'.")
    }
    if (Verbose) {
      cat("\nUsing user supplied decoy alignments.\n")
    }
    InputDecoys <- cbind(StaticDecoys,
                         "Response" = rep(0, nrow(StaticDecoys)))
    curr_sample <- ceiling(nrow(InputPairs) * DecoyScalar)
    if (curr_sample > nrow(StaticDecoys)) {
      curr_sample <- seq(nrow(StaticDecoys))
    } else {
      curr_sample <- sample(x = nrow(StaticDecoys),
                            size = curr_sample,
                            replace = FALSE)
    }
    InputDecoys <- InputDecoys[curr_sample, ]
    rownames(InputDecoys) <- NULL
    InputPairs <- cbind(InputPairs,
                        "Response" = rep(1, nrow(InputPairs)))
    EvalData <- rbind(InputPairs,
                      InputDecoys)
    
  } else {
    # neither, what should we do in this case?
    if (Verbose) {
      cat("\nNo decoy aligments supplied.\n")
    }
    EvalData <- cbind(InputPairs,
                      "Response" = rep(1, nrow(InputPairs)))
  }
  # if we got here, we should have a single `PairSummaries` object named EvalData
  # that has an extra column named 'Response'
  class(EvalData) <- c("data.frame",
                       "PairSummaries")
  attr(x = EvalData,
       which = "DefaultTranslationTable") <- attr(x = InputPairs,
                                                  which = "DefaultTranslationTable")
  attr(x = EvalData,
       which = "GeneCalls") <- attr(x = InputPairs,
                                    which = "GeneCalls")
  attr(x = EvalData,
       which = "KmerSize") <- attr(x = InputPairs,
                                   which = "KmerSize")
  
  # return(EvalData)
  candidate_data <- data.frame("Consensus" = EvalData$Consensus,
                               "FeatureDiff" = abs(EvalData$p1featurelength - EvalData$p2featurelength) / pmax(EvalData$p1featurelength,
                                                                                                               EvalData$p2featurelength),
                               "KDist" = EvalData$KDist,
                               "MaxMatch" = (EvalData$MaxMatch * 2L) / (EvalData$p1featurelength + EvalData$p2featurelength),
                               "Local_PID" = EvalData$Local_PID,
                               "Approx_Global_PID" = EvalData$Approx_Global_PID,
                               "Delta_Background" = EvalData$Delta_Background,
                               "Response" = EvalData$Response,
                               "Alignment" = EvalData$Alignment)
  
  if (EvaluationMethod == "glm") {
    # return a modeled response value
    int_res <- append_by_glm(Candidates = candidate_data)
    # return(list(EvalData,
    #             int_res))
    EvalData$criteria_value <- int_res$fitted.values
  } else if (EvaluationMethod == "lm") {
    # return a modeled value
    int_res <- append_by_lm(Candidates = candidate_data)
    # return(list(EvalData,
    #             int_res))
    EvalData$criteria_value <- int_res$fitted.values
  } else if (EvaluationMethod == "kmeans") {
    # return clusters
    int_res <- append_by_k(Candidates = candidate_data)
    EvalData$criteria_value <- int_res$cluster
    attr(x = EvalData,
         which = "centers") <- int_res$centers
    # return(EvalData)
  } else if (EvaluationMethod == "none") {
    # just return EvalData
    # return(EvalData)
  } else {
    stop("'EvaluationMethod' implies a method that is not supported.")
  }
  
  if (!is.null(FDRCriteria)) {
    if (!(names(FDRCriteria) %in% colnames(EvalData))) {
      stop("If 'FDRCriteria' is not NULL, it must be a named vector whose name corresponds to a column name for the supplied data.")
    }
    if (length(unique(EvalData$Response)) == 1L) {
      stop("'FDRCriteria' requires 'Response' data imparted by static, or generated decoy alignments.")
    }
    fdr_name <- names(FDRCriteria)[1L]
    fdr_lim <- unname(FDRCriteria)[1L]
    # true where is a decoy, because we're summing *up* the False Discovery Rate
    w1_logical <- EvalData$Response == 0
    w1_integer <- which(w1_logical)
    ranking <- order(EvalData[[fdr_name]],
                     decreasing = TRUE)
    rate <- cumsum(w1_logical[ranking]) / seq(nrow(EvalData))
    # the place in the ranking where we cross the fdr threshold
    w2_point <- min(which(rate >= fdr_lim))
    # drop every decoy alignment above the rate, and every candidate alignment
    # below the decoy rate limit
    retained <- ranking[seq(nrow(EvalData)) < w2_point]
    retained <- retained[!(retained %in% w1_integer)]
    
    EvalData <- EvalData[retained, ]
    rownames(EvalData) <- NULL
  }
  if (Verbose) {
    cat("\nCompleted!\n")
    tend <- Sys.time()
    print(tend - tstart)
  }
  return(EvalData)
}


