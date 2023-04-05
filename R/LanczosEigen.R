####
# Method to compute first k eigenvalues using Lanczos method
# Note that if k values are desired, m = 3k
# Originally written by Erik Wright
# 
# This method is not user-exposed
# Called in SuperTree() via dineof()
####
LanczosEigen <- function(A, m, tol=sqrt(.Machine$double.eps)) {
	if (m <= 1L)
		stop("m must be at least 2.")
	d <- nrow(A)
	if (m > d)
		stop("m can be at most ", d, ".")
	if (!isSymmetric(A))
		stop("A must be symmetric.")
	
	V <- matrix(NA_real_, d, m)
	V[, 1L] <- runif(d)
	
	V[, 1L] <- V[, 1L]/sqrt(sum(V[, 1L]^2))
	w <- A %*% V[, 1L]
	a <- t(w) %*% V[, 1L]
	w <- w - t(a %*% V[, 1L])
	for (j in 2:m) {
		B <- sqrt(sum(w^2))
		if (B < tol) { # avoid divide by ~zero
			V[, j] <- runif(d)
			V[, j] <- V[, j]/sqrt(sum(V[, j]^2))
		} else {
			V[, j] <- w/B
		}
		
		# enforce orthogonality
		for (i in seq_len(j - 1L)) {
			e <- sum(V[, j]*V[, i])
			if (abs(e) > tol)
				V[, j] <- V[, j] - e*V[, i]
		}
		V[, j] <- V[, j]/sqrt(sum(V[, j]^2))
		
		w <- A %*% V[, j]
		a <- t(w) %*% V[, j]
		w <- w - t(a %*% V[, j]) - B*V[, j - 1L]
	}
	
	TRI <- t(V) %*% A %*% V # tridiagonal
	# zero out off-tridiagonal, some entries stay due to inaccuracy
	TRI[abs(row(TRI) - col(TRI)) > 1] <- 0
	#return(TRI)
	e <- eigen(TRI, symmetric=TRUE)
	e$vectors <- V %*% e$vectors # eigenvectors of A
	o <- order(e$values, decreasing=TRUE)
	e$values <- e$values[o]
	e$vectors <- e$vectors[, o]
	e
}
