PlotNetworkInteger <- function(AdjMatrix, useRows=TRUE, useLabs=TRUE,
                        estimateEigens=(nrow(AdjMatrix)>24L),
                        pch=46, ...){
  diag(AdjMatrix) <- 0L
  if(useRows)
    rs <- rowSums(AdjMatrix)
  else
    rs <- colSums(AdjMatrix)
  lineplot <- which(AdjMatrix > 0 & lower.tri(AdjMatrix), arr.ind=TRUE)
  AdjMatrix <- -1*AdjMatrix
  diag(AdjMatrix) <- rs
  n <- nrow(AdjMatrix)
  if(estimateEigens){
    eigens <- LanczosEigen(AdjMatrix, 9L)
    vecs <- eigens$vectors[,c(7L,8L)]
  } else {
    eigens <- eigen(AdjMatrix, symmetric=TRUE)
    vecs <- eigens$vectors[,c(n-2L, n-1L)]
  }
  
  xv <- vecs[,1]
  yv <- vecs[,2]
  xv <- xv - min(xv)
  xv <- xv / max(xv)
  yv <- yv - min(yv)
  yv <- yv / max(yv)
  #print(lineplot)
  plot(NULL, xlab='Eigen1', ylab='Eigen2', 
       xaxt='n', yaxt='n', xlim=c(-0.05,1.05), ylim=c(-0.05,1.05), ...)
  for(i in seq_len(nrow(lineplot))){
    pos <- lineplot[i,]
    lines(x=xv[pos], y=yv[pos], col='gray', ...)
  }
  if(useLabs){
    if(is.null(rownames(AdjMatrix))){
      rn <- as.character(seq_len(n))
    } else {
      rn <- rownames(AdjMatrix)
    }
    text(xv, yv, labels=rn)
  } else {
    points(x=xv, y=yv, pch=pch, ...)
  }
  invisible(lineplot)
}

PlotContinuousNetwork <- function(Scores, Cutoff, 
                                  useRows=TRUE, useLabs=TRUE,
                                  estimateEigens=(nrow(AdjMatrix)>24L),
                                  largerEigens=FALSE, pch=46, 
                                  cols=colorRampPalette(c('#D81B60','white','#1E88E5'))(100),
                                  ...){
  Scores[] <- Scores - min(Scores, na.rm=TRUE)
  Scores[] <- Scores / max(Scores, na.rm=TRUE)
  Scores[is.na(Scores)] <- 0
  AdjMatrix <- Scores
  AdjMatrix[] <- as.integer(AdjMatrix >= Cutoff)
  diag(AdjMatrix) <- 0L
  if(useRows)
    rs <- rowSums(AdjMatrix)
  else
    rs <- colSums(AdjMatrix)
  lineplot <- which(AdjMatrix > 0 & lower.tri(AdjMatrix), arr.ind=TRUE)
  AdjMatrix <- -1*AdjMatrix
  diag(AdjMatrix) <- rs
  n <- nrow(AdjMatrix)
  if(estimateEigens){
    eigens <- LanczosEigen(AdjMatrix, 9L)
    #vecs <- eigens$vectors[,c(7L,8L)]
    vecs <- eigens$vectors[,eigens$values>1e-8]
  } else {
    eigens <- eigen(AdjMatrix, symmetric=TRUE)
    vecs <- eigens$vectors[,eigens$values>1e-8]
  }
  if(largerEigens)
    vecs <- vecs[,c(1,2)]
  else
    vecs <- vecs[,c(ncol(vecs)-1, ncol(vecs))]
  xv <- vecs[,1]
  yv <- vecs[,2]
  xv <- xv - min(xv)
  xv <- xv / max(xv)
  yv <- yv - min(yv)
  yv <- yv / max(yv)
  #print(lineplot)
  plot(NULL, xlab='Eigen1', ylab='Eigen2', 
       xaxt='n', yaxt='n', xlim=c(-0.05,1.05), ylim=c(-0.05,1.05), ...)
  
  for(i in seq_len(nrow(lineplot))){
    pos <- lineplot[i,]
    w <- Scores[pos[1], pos[2]]
    w <- (w - Cutoff) / (1-Cutoff)
    w <- round(w, 2)
    lines(x=xv[pos], y=yv[pos], col=cols[w*100], lwd=w, ...)
  }
  if(useLabs){
    if(is.null(rownames(AdjMatrix))){
      rn <- as.character(seq_len(n))
    } else {
      rn <- rownames(AdjMatrix)
    }
    text(xv, yv, labels=rn)
  } else {
    points(x=xv, y=yv, pch=pch, ...)
  }
}
