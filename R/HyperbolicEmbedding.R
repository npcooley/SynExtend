# Hyperbolic embedding using the Hydra algorithm from 
# https://arxiv.org/pdf/1903.08977.pdf
HyperbolicEmbedding <- function(M, d=2L, k=2, lambda=0.5){
  if(lambda < 0 || lambda > 1)
    stop("'lambda' must be between 0 and 1")
  if(!is.matrix(M))
    stop("'M' must be a matrix")
  if(!is.numeric(M))
    stop("'M' must be numeric or integer")
  if(!isSymmetric(M))
    stop("'M' must be a symmetric matrix!")
  A <- cosh(sqrt(k) * M)
  n <- nrow(A)
  # eigenvalues are greatest to smallest
  if(d < k / 6L){
    eigens <- LanczosEigen(A, 3L*d)
  } else {
    # no reason to use Lanczos if approximating more 
    # than half of the total eigenvalues
    eigens <- eigen(A, symmetric = TRUE)
    o <- order(eigens$values, decreasing=TRUE)
    eigens$values <- eigens$values[o]
    eigens$vectors <- eigens$vectors[, o]
  }
  
  # Create hyperbolic embedding
  X <- matrix(NA_real_, nrow=n, ncol=d+1L)
  colnums <- c(1L, seq(3*d-d,3*d))
  for(i in seq_len(d+1L)){
    if(i == 1){
      v <- sqrt(eigens$values[i]) * eigens$vectors[,i]
    } else {
      v <- sqrt(max(-1*eigens$values[colnums[i]], 0)) * eigens$vectors[,colnums[i]]
    }
    X[,i] <- v
  }
  
  # Project onto poincare ball
  u <- matrix(NA_real_, nrow=n, ncol=d)
  r <- numeric(n)
  s <- seq(2L, d+1)
  xm <- min(1L, X[,1L])
  for(i in seq_len(n)){
    u[i,] <- (X[i,s]) / sqrt(sum((X[i,s])**2))
    r[i] <- sqrt((X[i,1L] - xm) / (X[i,1L] + xm))
  }
  
  if(lambda > 0){
    const <- 2*pi / n
    thetas <- atan2(u[,2],u[,1])
    ranks <- rank(thetas, ties.method = 'first')
    thetas <- lambda*thetas + (1-lambda)*(ranks-1)*const
    u[,1L] <- r*cos(thetas)
    u[,2L] <- r*sin(thetas)
  } else {
    u <- u*r
  }
  
  rownames(u) <- rownames(M)
  colnames(u) <- paste0('V', seq_len(d))
  return(u)
}
