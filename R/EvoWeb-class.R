##### -- EvoWeb Class for output of EvoWeaver predictions --------------------
# author: Aidan Lakshman
# contact: ahl27@pitt.edu


## Defining S3 class for EvoWeb Object
##  This is just an S3 class to make output of EvoWeaver predictions 
##  easier to manage.

plot.EvoWeb <- function(x, NumSims=10, Gravity=0.05, Coulomb=0.1, Connection=5,
                         MoveRate=0.25, Cutoff=0.2, 
                         ColorPalette=topo.colors, Verbose=TRUE, ...){
  ## Springy embedding
  #
  # Uses 3 forces:
  # - Gravity (force towards (0,0))
  # - Coulomb force (repels nodes)
  # - Connection force
  web <- as.matrix(x)
  starttime <- Sys.time()
  if (Verbose) cat('Finding a descriptive embedding...\n')
  # Change similarity scores to distances
  web <- abs(web)
  web <- 1 - web
  web[web<Cutoff] <- 0
  
  # initialize to random values
  embedding <- matrix(runif(2*nrow(web)), ncol=2, nrow=nrow(web))
  true_dists <- cur_dists <- unclass(web)
  if (Verbose) pb <- txtProgressBar(max=NumSims, style=3)
  for (iter in seq_len(NumSims)){
    # Update distances
    for ( i in seq_len(nrow(cur_dists)) )
      for ( j in i:nrow(cur_dists) )
        cur_dists[i,j] <- cur_dists[j,i] <- sqrt(sum((embedding[i,]-embedding[j,])**2))
    
    # Gravity force
    g_vec <- -1 * Gravity * embedding
    
    # Coulomb and Attractive
    node_node <- g_vec
    node_node[] <- 0
    for ( i in seq_len(nrow(embedding)) ){
      for ( j in seq_len(nrow(embedding)) ){
        dir_force <- embedding[j,] - embedding[i,]
        r <- ifelse(cur_dists[i,j]==0, 1, cur_dists[i,j]**2)
        attr <- dir_force * (cur_dists[i,j] - true_dists[i,j]) * Connection
        repel <- -1 * Coulomb * dir_force / r
        node_node[i,] <- node_node[i,] + (attr + repel) / nrow(embedding)
      }
    }
    embedding <- embedding + MoveRate*(node_node + g_vec)
    if (Verbose) setTxtProgressBar(pb, iter)
  }
  
  center <- c(mean(embedding[,1]), mean(embedding[,2]))
  colsvec <- vapply(seq_len(nrow(embedding)), 
                    function(y) sqrt(sum((embedding[y,] - center)**2)),
                    FUN.VALUE=numeric(1))
  colors <- ColorPalette(length(colsvec))
  embedding <- embedding[order(colsvec, decreasing=FALSE),]
  if (Verbose)
    cat('\nDone.\n\nTime difference of', 
        round(difftime(Sys.time(), starttime, units = 'secs'), 2),
        'seconds.\n')
  plot(embedding[,1], embedding[,2], pch=19, cex=0.5, 
       xaxt='n', yaxt='n', ylab='', xlab='', col=colors,
       main='Force-directed embedding of COGs', ...)
}

summary.EvoWeb <- function(object, ...){
  cat('a EvoWeb object.\n')
  a <- attributes(object)
  cat('\tMethod used:', a$method, '\n')
  numGenes <- length(a$NAMES)
  numPreds <- sum(!is.na(object))
  cat('\tNumber of genes:', numGenes, '\n')
  cat('\tNumber of predictions:', numPreds, '\n')
  if ('model' %in% names(a)){
    cat('\tEnsemble model summary:\n')
    summary(a$model)
  }
}

show.EvoWeb <- function(object){
  cat('a EvoWeb object.\n')
  a <- attributes(object)
  cat('\tMethod used:', a$method, '\n')
  numGenes <- length(a$NAMES)
  numPreds <- sum(!is.na(object))
  cat('\tNumber of genes:', numGenes, '\n')
  cat('\tNumber of predictions:', numPreds, '\n')
  cat('\tPredictions:', numPreds, '\n\n')
  NextMethod()
}

print.EvoWeb <- function(x, ...){
  summary(x)
  cat('\nPredictions:\n\n')
  NextMethod()
}


as.data.frame.EvoWeb <- function(x, ...) {
  # Fall through to 'simMat' function
  NextMethod()
}

as.matrix.EvoWeb <- function(x, ...) {
  # Fall through to 'simMat' function
  NextMethod()
}

EvoWebMatToDf <- function(preds){
  stopifnot(is(preds,'matrix'))
  pair_locs <- upper.tri(preds)
  pairnames <- which(pair_locs, arr.ind=TRUE)
  pairentry1 <- rownames(preds)[pairnames[,'row']]
  pairentry2 <- colnames(preds)[pairnames[,'col']]
  AdjDf <- data.frame(Gene1=pairentry1, Gene2=pairentry2)
  AdjDf[,'Prediction'] <- preds[pair_locs]
  
  nc <- ncol(AdjDf)
  rtk <- vapply(seq_len(nrow(AdjDf)), function(x) sum(is.na(AdjDf[x,3:nc])) < (nc/2 - 1),
                FUN.VALUE=TRUE)
  AdjDf <- AdjDf[rtk,]
  rownames(AdjDf) <- NULL
  return(AdjDf)
}
