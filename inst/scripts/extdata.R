###### -- External data for SynExtend  ----------------------------------------
# Author: Nicholas Cooley
# email: npc19@pitt.edu
# this pseudocode relies on the ncbi command line utilities, they can be found
# here: https://www.ncbi.nlm.nih.gov/books/NBK179288/
# they must be installed, and R must have access to the executables
# this data was generated on 2022 04 21
# Self note:
# UPDATE DATALIST MANUALLY
# data was regenerated on 2022 09 22

suppressMessages(library(SynExtend))

TODAYSDATE <- paste0(unlist(strsplit(x = as.character(Sys.time()),
                                     split = "-| ")[[1]][1:3]),
                     collapse = "")

###### -- Entrez --------------------------------------------------------------

EntrezQuery <- paste("esearch -db assembly ",
                     "-query '",
                     "endosymbiont[All Fields] ",
                     'AND "complete genome"[filter] ',
                     'AND "RefSeq has annotation"[properties] ',
                     "NOT anomalous[filter]' ",
                     '| ',
                     'esummary ',
                     '| ',
                     'xtract -pattern DocumentSummary -element FtpPath_RefSeq',
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

###### -- Data import ---------------------------------------------------------
# save off GFFs as external non-R data
# save off `GeneCalls` as an object for examples

# save off one GFF for `gffToDataFrameExample`
CURLCOMMAND <- paste0("curl --output ",
                      paste0("~/Packages/SynExtend/inst/extdata/",
                             unlist(regmatches(x = GFFs[1],
                                               m = gregexpr(pattern = "[^/]+\\.gff\\.gz",
                                                            text = GFFs[1])))),
                      " ",
                      GFFs[1])

system(command = CURLCOMMAND,
       intern = FALSE)

Endosymbionts_GeneCalls <- vector(mode = "list",
                                    length = length(GFFs))

VignetteDB <- "~/Packages/SynExtend/inst/extdata/Endosymbionts.sqlite"

for (m1 in seq_along(GFFs)) {
  Endosymbionts_GeneCalls[[m1]] <- gffToDataFrame(GFF = GFFs[m1],
                                                    Verbose = TRUE)
  Seqs2DB(seqs = FNAs[m1],
          type = "FASTA",
          dbFile = VignetteDB,
          identifier = as.character(m1),
          verbose = TRUE)
}

names(Endosymbionts_GeneCalls) <- seq(length(Endosymbionts_GeneCalls))

Endosymbionts_Synteny <- FindSynteny(dbFile = VignetteDB,
                                       verbose = TRUE)

save(Endosymbionts_Synteny,
     file = "~/Packages/SynExtend/data/Endosymbionts_Synteny.RData",
     compress = "xz")

save(Endosymbionts_GeneCalls,
     file = "~/Packages/SynExtend/data/Endosymbionts_GeneCalls.RData",
     compress = "xz")

###### -- NucleotideOverlap ---------------------------------------------------

Endosymbionts_LinkedFeatures <- NucleotideOverlap(SyntenyObject = Endosymbionts_Synteny,
                                                  GeneCalls = Endosymbionts_GeneCalls,
                                                  Verbose = TRUE)

save(Endosymbionts_LinkedFeatures,
     file = "~/Packages/SynExtend/data/Endosymbionts_LinkedFeatures.RData",
     compress = "xz")

###### -- PairSummaries -------------------------------------------------------

Endosymbionts_Pairs01 <- PairSummaries(SyntenyLinks = Endosymbionts_LinkedFeatures,
                                         PIDs = TRUE,
                                         Score = TRUE,
                                         DBPATH = VignetteDB,
                                         Verbose = TRUE)

save(Endosymbionts_Pairs01,
     file = "~/Packages/SynExtend/data/Endosymbionts_Pairs01.RData",
     compress = "xz")

###### -- BlockExpansion ------------------------------------------------------

Endosymbionts_Pairs02 <- BlockExpansion(Pairs = Endosymbionts_Pairs01,
                                        NewPairsOnly = FALSE,
                                        DBPATH = VignetteDB,
                                        Verbose = TRUE)

save(Endosymbionts_Pairs02,
     file = "~/Packages/SynExtend/data/Endosymbionts_Pairs02.RData",
     compress = "xz")

###### -- BlockReconciliation -------------------------------------------------

Endosymbionts_Pairs03 <- BlockReconciliation(Pairs = Endosymbionts_Pairs02,
                                             ConservativeRejection = FALSE,
                                             Verbose = TRUE)
save(Endosymbionts_Pairs03,
     file = "~/Packages/SynExtend/data/Endosymbionts_Pairs03.RData",
     compress = "xz")

###### -- DisjointSet ---------------------------------------------------------

Endosymbionts_Sets <- DisjointSet(Pairs = Endosymbionts_Pairs03,
                                  Verbose = TRUE)

save(Endosymbionts_Sets,
     file = "~/Packages/SynExtend/data/Endosymbionts_Sets.RData",
     compress = "xz")

###### -- ExtractBy ----------------------------------------------------------- 
# no functions in the pipeline beyond this function
# no need to save off this for examples
# Endosymbionts_Gene_Communities <- ExtractBy(x = Endosymbionts_Pairs03,
#                                             y = VignetteDB,
#                                             z = Endosymbionts_Sets,
#                                             Verbose = TRUE)
# 
# save(Endosymbionts_Gene_Communities,
#      file = "~/Packages/SynExtend/data/Endosymbionts_Gene_Communities.RData",
#      compress = "xz")

