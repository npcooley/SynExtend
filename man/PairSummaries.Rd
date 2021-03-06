\name{PairSummaries}
\alias{PairSummaries}
\title{
Summarize connected pairs in a LinkedPairs object
}
\description{
Takes in a ``LinkedPairs'' object and gene calls, and returns a data.frame of paired features.
}
\usage{
PairSummaries(SyntenyLinks,
              DBPATH,
              PIDs = FALSE,
              IgnoreDefaultStringSet = FALSE,
              Verbose = FALSE,
              Model = "Generic",
              DefaultTranslationTable = "11",
              AcceptContigNames = TRUE,
              OffSetsAllowed = 2L,
              Storage = 1,
              ...)
}
\arguments{
  \item{SyntenyLinks}{
A \code{LinkedPairs} object. In previous versions of this function, a \code{GeneCalls} object was also required, but this object is now carried forward from \code{NucleotideOverlap} inside the \code{LinkedPairs} object.
}
  \item{DBPATH}{
A SQLite connection object or a character string specifying the path to the database file. Constructed from DECIPHER's \code{Seqs2DB} function. This path is always required as ``PairsSummaries'' computes the tetramer distance between paired sequences.
}
  \item{PIDs}{
Logical indicating whether to perform pairwise alignments. If \code{TRUE} all pairs will be aligned using DECIPHER's \code{AlignProfiles}. This step can be time consuming, especially for large numbers of pairs. Default is \code{FALSE}.
}
  \item{IgnoreDefaultStringSet}{
Logical indicating alignment type preferences. If \code{FALSE} (the default) pairs that can be aligned in amino acid space will be aligned as an \code{AAStringSet}. If \code{TRUE} all pairs will be aligned in nucleotide space. For \code{PairSummaries} to align the translation of a pair of sequences, both sequences must be tagged as coding in the ``GeneCalls'' object, and be the correct width for translation.
}
  \item{Verbose}{
Logical indicating whether or not to display a progress bar and print the time difference upon completion.
}
  \item{...}{
Arguments to be passed to \code{AlignProfiles}, and \code{DistanceMatrix}.
}
  \item{Model}{
A character string specifying a model to use to predict PIDs without performing an alignment. By default this argument is ``Generic'' specifying a generic PID prediction model based on PIDs computed from a randomly selected set of genomes. Currently no other models are included. Users may also supply their own model of type ``glm'' if they so desire in the form of an RData file. This model will need to take in some, or of the columns of statistics per pair that PairSummaries supplies.
}
  \item{DefaultTranslationTable}{
A character used to set the default translation table for \code{translate}. Is passed to \code{getGeneticCode}. Used when no translation table is specified in the ``GeneCalls'' object.
}
  \item{AcceptContigNames}{
Match names of contigs between gene calls object and synteny object. Where relevant, the first white space and everything following are removed from contig names. If \code{TRUE}, PairSummaries assumes that the contigs at each position in the synteny object and ``GeneCalls'' object are in the same order. Is automatically set to \code{TRUE} when ``GeneCalls'' are of class ``GRanges''. Is currently \code{TRUE} by default.
}
  \item{OffSetsAllowed}{
Integer vector defaulting to ``2L'' that sets the gap size that is allowed to be filled, if gaps are queried. The default value queries gaps of size 1. If set to ``NULL'' no gaps are queried. Setting to ``c(2L, 3L)'' would query all gaps of size 1 and 2.
}
  \item{Storage}{
Numeric indicating the approximate size a user wishes to allow for holding \code{StringSet}s in memory to extract gene sequences, in ``Gigabytes''. The lower \code{Storage} is set, the more likely that \code{PairSummaries} will need to reaccess \code{StringSet}s when extracting gene sequences. The higher \code{Storage} is set, the more sequences \code{PairSummaries} will attempt to hold in memory, avoiding the need to re-access the source database many times. Set to 1 by default, indicating that \code{PairSummaries} can store a ``Gigabyte'' of sequences in memory at a time.
}
}
\details{
The \code{LinkedPairs} object generated by \code{NucleotideOverlap} is a container for raw data that describes possible orthologous relationships, however ultimate assignment of orthology is up to user discretion. \code{PairSummaries} generates a clear table with relevant statistics for a user to work with as they choose. The option to align all pairs, though onerous can allow users to apply a hard threshold to predictions by PID, while built in models can allow more expedient thresholding from predicted PIDs.
}
\value{
A data.frame of class ``data.frame'' and ``PairSummaries'' of paired genes that are connected by syntenic hits. Contains columns describing the k-mers that link the pair. Columns ``p1'' and ``p2'' give the location ids of the the genes in the pair in the form ``DatabaseIdentifier_ContigIdentifier_GeneIdentifier''. ``ExactMatch'' provides an integer representing the exact number of nucleotides contained in the linking k-mers. ``TotalKmers'' provides an integer describing the number of distinct k-mers linking the pair. ``MaxKmer'' provides an integer describing the largest k-mer that links the pair. A column titled ``Concensus'' provides a value between zero and 1 indicating whether the kmers that link a pair of features are in the same position in each feature, with 1 indicating they are in exactly the same position and 0 indicating they are in as different a position as is possible. The ``Adjacent'' column provides an integer value ranging between 0 and 2 denoting whether a feature pair's direct neighbors are also paired. Gap filled pairs neither have neighbors, or are included as neighbors. The ``TetDist'' column provides the euclidean distance between oligonucleotide - of size 4 - frequences between predicted pairs. ``PIDType'' provides a character vector with values of ``NT'' where either of the pair indicates it is not a translatable sequence or ``AA'' where both sequences are translatable. If users choose to perform pairwise alignments there will be a ``PID'' column providing a numeric describing the percent identity between the two sequences. If users choose to predict PIDs using their own, or a provided model, a ``PredictedPID'' column will be provided.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{FindSynteny}}, \code{\link{Synteny-class}}
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
}
}
