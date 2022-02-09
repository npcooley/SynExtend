###### -- External data for SynExtend  ----------------------------------------
# Author: Nicholas Cooley
# email: npc19@pitt.edu
# this pseudocode relies on the ncbi command line utilities, they can be found
# here: https://www.ncbi.nlm.nih.gov/books/NBK179288/
# they must be installed, and R must have access to the executables
# this data was generated on 2022 02 09
# Self note:
# UPDATE DATALIST MANUALLY

suppressMessages(library(SynExtend))

TODAYSDATE <- paste0(unlist(strsplit(x = as.character(Sys.time()),
                                     split = "-| ")[[1]][1:3]),
                     collapse = "")

###### -- Entrez --------------------------------------------------------------

EntrezQuery <- paste("esearch -db assembly ",
                     "-query '",
                     "Nitrosocosmicus[organism] ",
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

Nitrosocosmicus_GeneCalls <- vector(mode = "list",
                                    length = length(GFFs))

VignetteDB <- "~/Packages/SynExtend/inst/extdata/Nitrosocosmicus.sqlite"

for (m1 in seq_along(GFFs)) {
  Nitrosocosmicus_GeneCalls[[m1]] <- gffToDataFrame(GFF = GFFs[m1],
                                                    Verbose = TRUE)
  Seqs2DB(seqs = FNAs[m1],
          type = "FASTA",
          dbFile = VignetteDB,
          identifier = as.character(m1),
          verbose = TRUE)
}

names(Nitrosocosmicus_GeneCalls) <- seq(length(Nitrosocosmicus_GeneCalls))

Nitrosocosmicus_Synteny <- FindSynteny(dbFile = VignetteDB,
                                       verbose = TRUE)

save(Nitrosocosmicus_Synteny,
     file = "~/Packages/SynExtend/data/Nitrosocosmicus_Synteny.RData",
     compress = "xz")

save(Nitrosocosmicus_GeneCalls,
     file = "~/Packages/SynExtend/data/Nitrosocosmicus_GeneCalls.RData",
     compress = "xz")

###### -- NucleotideOverlap ---------------------------------------------------

Nitrosocosmicus_LinkedFeatures <- NucleotideOverlap(SyntenyObject = Nitrosocosmicus_Synteny,
                                                    GeneCalls = Nitrosocosmicus_GeneCalls,
                                                    Verbose = TRUE)

save(Nitrosocosmicus_LinkedFeatures,
     file = "~/Packages/SynExtend/data/Nitrosocosmicus_LinkedFeatures.RData",
     compress = "xz")

###### -- PairSummaries -------------------------------------------------------

Nitrosocosmicus_Pairs01 <- PairSummaries(SyntenyLinks = Nitrosocosmicus_LinkedFeatures,
                                         PIDs = TRUE,
                                         DBPATH = VignetteDB,
                                         Verbose = TRUE)

save(Nitrosocosmicus_Pairs01,
     file = "~/Packages/SynExtend/data/Nitrosocosmicus_Pairs01.RData",
     compress = "xz")

###### -- BlockExpansion ------------------------------------------------------

Nitrosocosmicus_Pairs02 <- BlockExpansion(Pairs = Nitrosocosmicus_Pairs01,
                                          NewPairsOnly = FALSE,
                                          DBPATH = VignetteDB,
                                          Verbose = TRUE)

save(Nitrosocosmicus_Pairs02,
     file = "~/Packages/SynExtend/data/Nitrosocosmicus_Pairs02.RData",
     compress = "xz")

###### -- BlockReconciliation -------------------------------------------------

Nitrosocosmicus_Pairs03 <- BlockReconciliation(Pairs = Nitrosocosmicus_Pairs02,
                                               ConservativeRejection = FALSE,
                                               Verbose = TRUE)
save(Nitrosocosmicus_Pairs03,
     file = "~/Packages/SynExtend/data/Nitrosocosmicus_Pairs03.RData",
     compress = "xz")

###### -- DisjointSet ---------------------------------------------------------

Nitrosocosmicus_Sets <- DisjointSet(Pairs = Nitrosocosmicus_Pairs03,
                                    Verbose = TRUE)

save(Nitrosocosmicus_Sets,
     file = "~/Packages/SynExtend/data/Nitrosocosmicus_Sets.RData",
     compress = "xz")

###### -- ExtractBy -----------------------------------------------------------

# no functions past extractby, no need to create package data from it

###### -- other stuff ---------------------------------------------------------

# some very arbitrary subsetting
w1 <- (Nitrosocosmicus_Pairs03$ExactMatch * 2L) / (Nitrosocosmicus_Pairs03$p1FeatureLength + Nitrosocosmicus_Pairs03$p2FeatureLength)
w2 <- w1 > 0.2 & w1 != 0
w3 <- Nitrosocosmicus_Pairs03$Consensus >= 0.6
w4 <- abs(Nitrosocosmicus_Pairs03$p1FeatureLength - Nitrosocosmicus_Pairs03$p2FeatureLength) / apply(X = Nitrosocosmicus_Pairs03[, c("p1FeatureLength", "p2FeatureLength")],
                                                                                                     MARGIN = 1L,
                                                                                                     FUN = max) < 0.5

Nitrosocosmicus_Pairs04 <- Nitrosocosmicus_Pairs03[w2 & w3 & w4, ]

suppressMessages(library(igraph))

w5 <- Nitrosocosmicus_Pairs04$p1 %in% Nitrosocosmicus_Sets[[which.max(lengths(Nitrosocosmicus_Sets))]] &
  Nitrosocosmicus_Pairs04$p2 %in% Nitrosocosmicus_Sets[[which.max(lengths(Nitrosocosmicus_Sets))]]
# what is the largest community
g <- graph_from_data_frame(d = Nitrosocosmicus_Pairs04[w5, c("p1", "p2", "PID")],
                           directed = FALSE)
# does it have distinct subcommunities
sg <- cluster_louvain(graph = g)
plot(sg,
     g,
     vertex.label = NA)

# extract and align and look at the dendrogram
x <- ExtractBy(x = Nitrosocosmicus_Pairs04,
               y = Nitrosocosmicus_Sets[which.max(lengths(Nitrosocosmicus_Sets))],
               Method = "clusters",
               Verbose = TRUE,
               DBPATH = VignetteDB)

ali <- AlignSeqs(myXStringSet = x[[1]])

IdClusters(myDistMatrix = DistanceMatrix(ali,
                                         includeTerminalGaps = TRUE),
           method = "NJ",
           type = "dendrogram",
           showPlot = TRUE)






