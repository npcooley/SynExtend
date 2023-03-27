####
# Method to impute missing values for distance matrices
# Uses Lanczos method for eigen values 
# Written by Aidan Lakshman
# 
# This method is not user exposed
# Called in SuperTree()
####
dineof <- function(X, start=1, 
                   num_iter = min(5, ncol(X)%/%10), tol=1e-4,
                   verbose = FALSE){
  if(verbose){
    cat("  Imputing missing distances...")
  }
  if (start <= 0)
    start <- 1
  
  num_iter <- start+num_iter-1
  n <- ncol(X) 
  
  # Flag missing values
  lowertrimult <- lower.tri(X)
  pos <- lowertrimult & is.na(X)
  posmissing <- which(pos)
  pospresent <- which(lowertrimult & !is.na(X))
  
  # Set aside some for cross validation
  cv_size <- min((0.01*n*n+40)/2, (0.03*n*n)/2)
  cv_size <- floor(cv_size)
  
  cv_set <- sample(pospresent, cv_size)
  cv_vals <- X[cv_set]
  X[cv_set] <- NA_real_
  
  all_to_update <- c(posmissing, cv_set)
  rev_atu <- all_to_update - 1
  rev_atu <- (rev_atu %% n)*n + (rev_atu %/% n) + 1
  
  # Nick Christoffersen says fill with zeros to ensure positive semidefinite
  # Random matrix eigenvalues tightly bounded around sqrt(n)
  X[is.na(X)] <- 0
  
  make_sym <- function(m){
    m <- m*lowertrimult
    t(t(m) + m) 
  }
  # Make matrix symmetric
  X <- X * lowertrimult
  X <- make_sym(X)
  best_eof <- n
  best_X <- cur_X <- X
  prev_e <- cur_e <- new_e <- best_e <- Inf
  eof <- start-1
  
  ## Make matrix symmetric
  
  if (verbose){
    progressstar <- c(" | "," / "," - "," \\ ")
    psl <- length(progressstar)
    
    waitinc <- 0.1
    catcounter <- 0
    cat("\n  Optimizing: ")
    firststring <- "\r  Optimizing:"
    lastprint <- Sys.time()
  }
  while(num_iter > eof || cur_e < prev_e){
    # Increment number of EOF, reset values
    eof <- eof + 1
    secondstring <- paste0("(Testing ", eof, ' of ', max(num_iter,eof), ' EOFs)')
    cur_X <- X
    cur_e <- 0
    new_e <- Inf
    prev_e <- cur_e
    
    while(abs(cur_e - new_e) > tol){
      if(verbose){
        if(Sys.time() - lastprint > waitinc){
          if(!is.infinite(new_e) && !is.infinite(cur_e)){
            cat(firststring, secondstring, progressstar[(catcounter)+1], "Loss:", 
                format(cur_e,digits=4,nsmall=4), 
                rep(' ', 20))
          } else {
            cat(firststring, secondstring, progressstar[(catcounter)+1],
                "Loss: ", rep(' ', 20))
          }
          catcounter <- (catcounter+1) %% psl
          lastprint <- Sys.time()
        }
      }
      
      # Repeat SVD until convergence
      cur_e <- new_e
      s <- LanczosEigen(cur_X, min(eof*3,n))
      eigenvalues <- s$values[seq_len(eof)]
      eigenvecs <- s$vectors[,seq_len(eof), drop=FALSE]
      
      preds <- eigenvecs %*% (t(eigenvecs) * eigenvalues)
      preds <- preds[all_to_update]
      cur_X[all_to_update] <- preds
      cur_X[rev_atu] <- preds
      new_e <- sqrt(mean((cur_X[cv_set] - cv_vals)**2))
      if(new_e > cur_e + 10*tol) break
    }
    cur_e <- new_e
    if(cur_e < best_e){
      best_e <- cur_e
      best_eof <- eof
      best_X <- cur_X
    }
  }
  
  if(verbose){
    cat("\r  Optimizing: Done.", rep(' ', 25), '\n')
  }
  
  # reset the cross validation set
  best_X[cv_set] <- cv_vals
  best_X <- make_sym(best_X)
  
  return(list(X=best_X, eof=best_eof, err=best_e))
}
