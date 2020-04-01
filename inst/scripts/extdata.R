###### -- External data for SynExtend  ----------------------------------------
# Author: Nicholas Cooley
# email: npc19@pitt.edu
# this pseudocode relies on the ncbi command line utilities, they can be found
# here: https://www.ncbi.nlm.nih.gov/books/NBK179288/
# they must be installed, and R must have access to the executables
# this data was generated on 2020 02 18
# save commands and paths have been commented and should be adjusted to fit
# the needs of a user


###### -- Libraries -----------------------------------------------------------

# library(devtools)
# install_github(repo = "npcooley/Heron")
library(Heron)

###### -- Arguments -----------------------------------------------------------

# no command line arguments

###### -- Data ----------------------------------------------------------------

# no external data

###### -- Functions -----------------------------------------------------------

# no adhoc functions present

###### -- Code Body -----------------------------------------------------------

EntrezQuery <- paste("esearch -db assembly ",
                     "-query '",
                     "thaumarchaeota[organism] ",
                     'AND "complete genome"[filter] ',
                     'AND "genbank has annotation"[properties] ',
                     "NOT anomalous[filter]' ",
                     '| ',
                     'esummary ',
                     '| ',
                     'xtract -pattern DocumentSummary -element FtpPath_GenBank',
                     sep = "")

FtPPaths <- system(command = EntrezQuery,
                   intern = TRUE,
                   timeout = 300L)

FNAs <- unname(sapply(FtPPaths,
                      function(x) paste(x,
                                        "/",
                                        strsplit(x,
                                                 split = "/",
                                                 fixed = TRUE)[[1]][10],
                                        "_genomic.fna.gz",
                                        sep = "")))

GFFs <- unname(sapply(FtPPaths,
                      function(x) paste(x,
                                        "/",
                                        strsplit(x,
                                                 split = "/",
                                                 fixed = TRUE)[[1]][10],
                                        "_genomic.gff.gz",
                                        sep = "")))


# NewDB <- "~/Dropbox/Packages/SynExtend/inst/extdata/VignetteSeqs.sqlite"

FNAs <- FNAs[c(1, 12, 13)]
GFFs <- GFFs[c(1, 12, 13)]

for (m1 in seq_along(FNAs)) {
  X <- readDNAStringSet(filepath = FNAs[m1])
  X <- X[order(width(X),
               decreasing = TRUE)]
  
  Seqs2DB(seqs = X,
          type = "XStringSet",
          dbFile = NewDB,
          identifier = as.character(m1),
          verbose = TRUE)
}

GFFsAsDFs <- vector(mode = "list",
                    length = length(GFFs)) -> GFFsAsGRanges

for (m1 in seq_along(GFFsAsDFs)) {
  GFFsAsDFs[[m1]] <- gffToDataFrame(GFF = GFFs[m1],
                                    Verbose = TRUE)
  GFFsAsGRanges[[m1]] <- rtracklayer::import(GFFs[m1])
}

###### -- Save Objects --------------------------------------------------------

# no objects saved...
