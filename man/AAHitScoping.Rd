\name{AAHitScoping}
\alias{AAHitScoping}
\title{
Adjust the scope of kmer hits between feature and genome space.
}
\description{
This function is designed to work internally to functions within \code{SynExtend} so it works on relatively simple atomic vectors and has little overhead checking.
}
\usage{
AAHitScoping(hitlist,
             fstrand1,
             fstart1,
             fstop1,
             fstrand2,
             fstart2,
             fstop2)
}
\arguments{
  \item{hitlist}{
  A list containing matrices produced by \code{\link{SearchIndex}}.
}
  \item{fstrand1}{
  An integer vector of 0s and 1s describing the strand of features.
}
  \item{fstart1}{
  Integer; a vector of left bounds of features.
}
  \item{fstop1}{
  Integer; a vector of right bounds of features.
}
  \item{fstrand2}{
  An integer vector of 0s and 1s describing the strand of features.
}
  \item{fstart2}{
  Integer; a vector of left bounds of features.
}
  \item{fstop2}{
  Integer; a vector of right bounds of features.
}
}
\details{
  \code{AAHitScoping} converts the hits returned by \code{\link{SearchIndex}} from feature-to-feature context genome-to-genome context.
}
\value{
  A list of matrices.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{NucleotideOverlap}}, \code{\link{SummarizePairs}}
}
\examples{
# a contrived example
x <- list(matrix(c(1L, 3L, 1L, 3L,
                   5L, 7L, 5L, 7L),
                 nrow = 4L,
                 ncol = 2L),
          matrix(c(2L, 4L, 2L, 4L),
                 nrow = 4L,
                 ncol = 1L))

# Feature 1 (query): forward strand, genome positions 100–199
# Feature 2 (query): reverse strand, genome positions 300–449
vals01 <- c(0L, 1L)   # 0 = forward, 1 = reverse
vals02 <- c(100L, 300L)
vals03 <- c(199L, 449L)

## Subject features mirror the same layout
vals04 <- c(0L, 1L)
vals05 <- c(500L, 700L)
vals06 <- c(599L, 849L)

result <- AAHitScoping(hitlist  = x,
                       fstrand1 = vals01,
                       fstart1 = vals02,
                       fstop1 = vals03,
                       fstrand2 = vals04,
                       fstart2 = vals05,
                       fstop2 = vals06)
}
