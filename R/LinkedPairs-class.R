`[.LinkedPairs` <- function(x, i, j, ...) {
  ans <- NextMethod("[", x)
  
  if (missing(j))
    return(ans)
  
  d <- dim(x)
  I <- seq_len(d[1])
  J <- seq_len(d[1])
  d <- dimnames(x)
  names(I) <- d[[1]]
  names(J) <- d[[1]]
  I <- I[i]
  J <- J[j]
  if (length(I) >= 2 &&
      length(I)==length(J) &&
      all(I==J) &&
      !any(duplicated(I)))
    class(ans) <- "LinkedPairs"
  return(ans)
}

print.LinkedPairs <- function(x,
                            quote = FALSE,
                            right = TRUE,
                            ...) {
  d <- dim(x)
  if (is.null(d)) {
    stop ("x must be a square object of class 'LinkedPairs'.")
  }
  m <- matrix("",
              nrow = d[1],
              ncol = d[2],
              dimnames = dimnames(x))
  
  for (m1 in seq_len(d[1])) {
    for (m2 in seq_len(d[2])) {
      if (m1 > m2) {
        ######
        # Lower Triangle
        ######
        k <- sum(x[m2, m1][[1]][, "ExactOverlap"])
        m[m1, m2] <- paste(k,
                           ifelse(k == 1,
                                  "Nucleotide",
                                  "Nucleotides"),
                           sep = " ")
      } else if (m1 < m2) {
        ######
        # Upper Triangle
        ######
        k <- nrow(x[m1, m2][[1]])
        m[m1, m2] <- paste(k,
                           ifelse(k == 1,
                                  "Pair",
                                  "Pairs"),
                           sep = " ")
      } else if (m1 == m2 &
                 !is.null(x[m1, m2][[1]])) {
        k <- nrow(x[m1, m2][[1]])
        m[m1, m2] <- paste(k,
                           ifelse(k == 1,
                                  "Gene",
                                  "Genes"),
                           sep = " ")
      } else if (m1 == m2 &
                 is.null(x[m1, m2][[1]])) {
        m[m1, m2] <- ""
      }
    }
  }
  print(m,
        quote = quote,
        right = right,
        ...)
  invisible(x)
}