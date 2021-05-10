\name{DisjointSet}
\alias{DisjointSet}
\title{
Return single linkage clusters from \code{PairSummaries} objects.
}
\description{
Takes in a \code{PairSummaries} object and return a list of identifiers organized into single linkage clusters.
}
\usage{
DisjointSet(Pairs,
            Verbose = FALSE)
}
\arguments{
  \item{Pairs}{
A \code{PairSummaries} object.
}
  \item{Verbose}{
Logical indicating whether to print progress bars and messages. Defaults to \code{FALSE}.
}
}
\details{
Takes in a \code{PairSummaries} object and return a list of identifiers organized into single linkage clusters.
}
\value{
Returns a list of character vectors representing IDs of sequence features, typically genes.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{FindSynteny}}, \code{\link{Synteny-class}}, \code{\link{PairSummaries}}, \code{\link{FindSets}}
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
