\name{ExtractBy}
\alias{ExtractBy}
\title{
Extract and organize \code{XStringSet}s of sequences represented in a \code{PairSummaries} object.
}
\description{
Takes in a \code{PairSummaries} object and an optional vector of cluster representatives. Return an \code{XStringSet} of the sequences present in the \code{PairSummaries}, or when cluster representatives are provided, a list of \code{XStringSet}s of the sequences that make up the provided clusters.
}
\usage{
ExtractBy(x,
          y = NULL,
          DBPATH,
          Method = "all",
          DefaultTranslationTable = "11",
          Translate = TRUE,
          Storage = 1,
          Verbose = FALSE)
}
\arguments{
  \item{x}{
A \code{PairSummaries} object.
}
  \item{y}{
An optional \code{list} containing the ids of sequences in the \code{PairSummaries} object.
}
  \item{DBPATH}{
A SQLite connection object or a character string specifying the path to the database file. Constructed from DECIPHER's \code{Seqs2DB} function.
}
  \item{Method}{
How to extract sequences from the \code{PairSummaries} object. Currently only the methods ``all'' and ``clusters'' are supported.
}
  \item{Translate}{
If \code{TRUE} return \code{AAStringSet}s where possible.
}
  \item{DefaultTranslationTable}{
Currently Not Implemented! When implemented will allow for designation of a specific translation table if one is not indicated in the \code{GeneCalls} attribute of the \code{PairSummaries} object.
}
  \item{Storage}{
Numeric indicating the approximate size a user wishes to allow for holding \code{StringSet}s in memory to extract gene sequences, in ``Gigabytes''. The lower \code{Storage} is set, the more likely that \code{ExtractBy} will need to reaccess \code{StringSet}s when extracting gene sequences. The higher \code{Storage} is set, the more sequences \code{ExtractBy} will attempt to hold in memory, avoiding the need to re-access the source database many times. Set to 1 by default, indicating that \code{ExtractBy} can store a ``Gigabyte'' of sequences in memory at a time.
}
  \item{Verbose}{
Logical indicating whether to print progress bars and messages. Defaults to \code{FALSE}.
}
}
\details{
Takes in a \code{PairSummaries} object and an optional vector of cluster representatives. Return an \code{XStringSet} of the sequences present in the \code{PairSummaries}, or when cluster representatives are provided, a list of \code{XStringSet}s of the sequences that make up the provided clusters.
}
\value{
Returns either a \code{XStringSet} or a list of \code{XStringSet}s.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{FindSynteny}}, \code{\link{Synteny-class}}, \code{\link{PairSummaries}}, \code{\link{DisjointSet}}
}
\examples{
DBPATH <- system.file("extdata",
                      "VignetteSeqs.sqlite",
                      package = "SynExtend")
Syn <- FindSynteny(dbFile = DBPATH)
GeneCalls <- vector(mode = "list",
                    length = ncol(Syn))

GeneCalls[[1L]] <- gffToDataFrame(GFF = system.file("extdata",
                                                    "GCA_006740685.1_ASM674068v1_genomic.gff.gz",
                                                    package = "SynExtend"),
                                  Verbose = TRUE)
GeneCalls[[2L]] <- gffToDataFrame(GFF = system.file("extdata",
                                                    "GCA_000956175.1_ASM95617v1_genomic.gff.gz",
                                                    package = "SynExtend"),
                                  Verbose = TRUE)
GeneCalls[[3L]] <- gffToDataFrame(GFF = system.file("extdata",
                                                    "GCA_000875775.1_ASM87577v1_genomic.gff.gz",
                                                    package = "SynExtend"),
                                  Verbose = TRUE)
names(GeneCalls) <- seq(length(GeneCalls))
Links <- NucleotideOverlap(SyntenyObject = Syn,
                           GeneCalls = GeneCalls,
                           LimitIndex = FALSE,
                           Verbose = TRUE)
PredictedPairs <- PairSummaries(SyntenyLinks = Links,
                                DBPATH = DBPATH,
                                PIDs = FALSE,
                                AcceptContigNames = TRUE,
                                Verbose = TRUE)
PresentSeqs <- ExtractBy(x = PredictedPairs,
                         Method = "all",
                         DBPATH = DBPATH,
                         Verbose = TRUE)
Clusters <- DisjointSet(Pairs = PredictedPairs,
                        Verbose = TRUE)
SeqsByClusters <- ExtractBy(x = PredictedPairs,
                            y = Clusters,
                            Method = "clusters",
                            DBPATH = DBPATH,
                            Verbose = TRUE)

# Alternatively the same seqs can be accessed from the NCBI FTP site
# And gene calls can be accessed with the rtracklayer
\dontrun{
DBPATH <- tempfile()
FNAs <- c("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/006/740/685/GCA_006740685.1_ASM674068v1/GCA_006740685.1_ASM674068v1_genomic.fna.gz",
          "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/956/175/GCA_000956175.1_ASM95617v1/GCA_000956175.1_ASM95617v1_genomic.fna.gz",
          "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/875/775/GCA_000875775.1_ASM87577v1/GCA_000875775.1_ASM87577v1_genomic.fna.gz")
for (m1 in seq_along(FNAs)) {
 X <- readDNAStringSet(filepath = FNAs[m1])
 X <- X[order(width(X),
              decreasing = TRUE)]
 
 Seqs2DB(seqs = X,
         type = "XStringSet",
         dbFile = DBPATH,
         identifier = as.character(m1),
         verbose = TRUE)
  }

GeneCalls <- vector(mode = "list",
                    length = ncol(Syn))
GeneCalls[[1L]] <- rtracklayer::import(system.file("extdata",
                                                   "GCA_006740685.1_ASM674068v1_genomic.gff.gz",
                                                   package = "SynExtend"))
GeneCalls[[2L]] <- rtracklayer::import(system.file("extdata",
                                                   "GCA_000956175.1_ASM95617v1_genomic.gff.gz",
                                                   package = "SynExtend"))
GeneCalls[[3L]] <- rtracklayer::import(system.file("extdata",
                                                   "GCA_000875775.1_ASM87577v1_genomic.gff.gz",
                                                   package = "SynExtend"))
names(GeneCalls) <- seq(length(GeneCalls))
Links <- NucleotideOverlap(SyntenyObject = Syn,
                           GeneCalls = GeneCalls,
                           LimitIndex = FALSE,
                           Verbose = TRUE)
PredictedPairs <- PairSummaries(SyntenyLinks = Links,
                                DBPATH = DBPATH,
                                PIDs = FALSE,
                                AcceptContigNames = TRUE,
                                Verbose = TRUE)
PresentSeqs <- ExtractBy(x = PredictedPairs,
                         Method = "all",
                         DBPATH = DBPATH,
                         Verbose = TRUE)
Clusters <- DisjointSet(Pairs = PredictedPairs,
                        Verbose = TRUE)
SeqsByClusters <- ExtractBy(x = PredictedPairs,
                            y = Clusters,
                            Method = "clusters",
                            DBPATH = DBPATH,
                            Verbose = TRUE)
}
}
