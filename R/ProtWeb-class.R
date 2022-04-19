##### -- ProtWeb Class for output of ProtWeaver predictions --------------------
# author: Aidan Lakshman
# contact: ahl27@pitt.edu


## Defining S3 class for ProtWeb Object

#### DEFINITION ####
####################

#### METHODS ####

#################

getData <- function(x, ...) UseMethod('getData')

plot.ProtWeb <- function(web, nsim=10, gravity=0.05, coulomb=0.1, connection=5,
                         move_rate=0.25, cutoff=0.2, 
                         verbose=T, colorpalette=topo.colors, ...){
  ## Springy embedding
  #
  # Uses 3 forces:
  # - gravity (force towards (0,0))
  # - Coulomb force (repels nodes)
  # - Connection force
  
  if (verbose) cat('Finding a descriptive embedding...\n')
  # Change similarity scores to distances
  web <- abs(web)
  web <- 1 - web
  web[web<cutoff] <- 0
  
  # initialize to random values
  embedding <- matrix(runif(2*nrow(web)), ncol=2, nrow=nrow(web))
  true_dists <- cur_dists <- unclass(web)
  if (verbose) pb <- txtProgressBar(max=nsim, style=3)
  for (iter in 1:nsim){
    # Update distances
    for ( i in 1:nrow(cur_dists) )
      for ( j in i:nrow(cur_dists) )
        cur_dists[i,j] <- cur_dists[j,i] <- sqrt(sum((embedding[i,]-embedding[j,])**2))
    
    # Gravity force
    g_vec <- -1 * gravity * embedding
    
    # Coulomb and Attractive
    node_node <- g_vec
    node_node[] <- 0
    for ( i in 1:nrow(embedding) ){
      for ( j in 1:nrow(embedding) ){
        dir_force <- embedding[j,] - embedding[i,]
        r <- ifelse(cur_dists[i,j]==0, 1, cur_dists[i,j]**2)
        attr <- dir_force * (cur_dists[i,j] - true_dists[i,j]) * connection
        repel <- -1 * coulomb * dir_force / r
        node_node[i,] <- node_node[i,] + (attr + repel) / nrow(embedding)
      }
    }
    embedding <- embedding + move_rate*(node_node + g_vec)
    if (verbose) setTxtProgressBar(pb, iter)
  }
  
  #return(embedding)
  center <- c(mean(embedding[,1]), mean(embedding[,2]))
  colsvec <- sapply(1:nrow(embedding), function(x) sqrt(sum((embedding[x,] - center)**2)))
  colors <- colorpalette(length(colsvec))
  embedding <- embedding[order(colsvec, decreasing=F),]
  plot(embedding[,1], embedding[,2], pch=19, cex=0.5, 
       xaxt='n', yaxt='n', ylab='', xlab='', col=colors,
       main='Force-directed embedding of COGs', ...)
}

summary.ProtWeb <- function(x, ...){
  cat('a ProtWeb object.\n')
  a <- attributes(x)
  cat('\tMethod used:', a$method, '\n')
  d <- getData(x)
  numGenes <- ncol(d)
  numPreds <- sum(upper.tri(d,diag=TRUE) & !is.na(d))
  cat('\tNumber of genes:', numGenes, '\n')
  cat('\tNumber of predictions:', numPreds, '\n')
  if ('model' %in% names(a)){
    cat('\tEnsemble model summary:\n')
    summary(a$model)
  }
}

show.ProtWeb <- function(x, ...){
  summary(x)
}

getData.ProtWeb <- function(x, asDf=F, ...){
  dims <- dim(x)
  rnames <- rownames(x)
  cnames <- colnames(x)
  attributes(x) <- NULL
  arr <- array(x, dim=dims)
  rownames(arr) <- rnames
  colnames(arr) <- cnames
  if (asDf){
    arr <- ProtWebMatToDf(arr)
  }
  return(arr)
}

print.ProtWeb <- function(x, ...){
  summary(x)
}

ProtWebMatToDf <- function(preds){
  stopifnot(is(preds,'matrix'))
  pair_locs <- upper.tri(preds)
  pairnames <- which(pair_locs, arr.ind=T)
  pairentry1 <- rownames(preds)[pairnames[,'row']]
  pairentry2 <- colnames(preds)[pairnames[,'col']]
  AdjDf <- data.frame(Gene1=pairentry1, Gene2=pairentry2)
  AdjDf[,'Prediction'] <- preds[pair_locs]
  
  nc <- ncol(AdjDf)
  rtk <- vapply(1:nrow(AdjDf), function(x) sum(is.na(AdjDf[x,3:nc])) < (nc/2 - 1),
                FUN.VALUE=TRUE)
  AdjDf <- AdjDf[rtk,]
  rownames(AdjDf) <- NULL
  return(AdjDf)
}
