###### -- External data for SynExtend  ----------------------------------------
# Author: Nicholas Cooley
# email: nicholas.cooley@ul.ie
# this script relies on the ncbi command line utilities, they can be found
# here: https://www.ncbi.nlm.nih.gov/books/NBK179288/
# they must be installed, and R must have access to the executables

suppressMessages(library(SynExtend))
suppressMessages(library(DBI))
# functions that have been changed or added will need to be sourced explicitly
source(file = "~/Repos/SynExtend/R/FrameDownward.R",
       echo = FALSE)
source(file = "~/Repos/SynExtend/R/SquaregffBy.R",
       echo = FALSE)
source(file = "~/Repos/SynExtend/R/FeaturesFromDF.R",
       echo = FALSE)
source(file = "~/Repos/SynExtend/R/CreateDecoys.R",
       echo = FALSE)
source(file = "~/Repos/SynExtend/R/SummarizePairs.R",
       echo = FALSE)
source(file = "~/Repos/SynExtend/R/EvaluatePairs.R",
       echo = FALSE)
# functions to remove:
# gffToDataFrame
# ClusterByWhatever
# ExpandDiagonal
# a whole bunch of others, basically anything from my set that relies on an
# old version of example data

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

# keep the example data small...
if (length(FtPPaths) > 4L) {
  # setting the seed here is a little superfluous, as the data size will almost certainly
  # change before this is run again, or before anyone comes back to run it themselves
  set.seed(1986)
  FtPPaths <- sample(FtPPaths,
                     size = 4,
                     replace = FALSE)
}

# sometimes the ftp paths don't play nicely, replace with https if necessary:
# this might be a result of being in the EU, or it could be something else
# ¯\_(ツ)_/¯
FtPPaths <- sub(pattern = "^ftp",
                replacement = "https",
                x = FtPPaths)

adds <- mapply(SIMPLIFY = TRUE,
               USE.NAMES = FALSE,
               FUN = function(x, y) {
                 paste0(x,
                        "/",
                        y[10],
                        c("_genomic.fna.gz",
                          "_genomic.gff.gz",
                          "_protein.faa.gz"))
               },
               x = FtPPaths,
               y = strsplit(x = FtPPaths,
                            split = "/",
                            fixed = TRUE))
fnas <- adds[1, , drop = TRUE]
gffs <- adds[2, , drop = TRUE]
amns <- adds[3, , drop = TRUE]


###### -- Data import ---------------------------------------------------------
# save off gffs as external non-R data
# save off `GeneCalls` as an object for examples

# save off one GFF for `SquaregffBy's example`
CURLCOMMAND <- paste0("curl --output ",
                      paste0("~/Repos/SynExtend/inst/extdata/",
                             unlist(regmatches(x = gffs[1],
                                               m = gregexpr(pattern = "[^/]+\\.gff\\.gz",
                                                            text = gffs[1])))),
                      " ",
                      gffs[1])

system(command = CURLCOMMAND,
       intern = FALSE)

genecalls <- vector(mode = "list",
                    length = length(gffs))

VignetteDB01 <- "~/Repos/SynExtend/inst/extdata/example_db.sqlite"
VignetteDB02 <- tempfile()

for (m1 in seq_along(gffs)) {
  x <- readDNAStringSet(filepath = fnas[m1])
  if (sum(width(x)) > 1e6) {
    next
  }
  
  tmp_obj <- rtracklayer::import(con = gffs[m1])
  genecalls[[m1]] <- SquaregffBy(gff_object = tmp_obj,
                                               verbose = TRUE)
  Seqs2DB(seqs = x,
          type = "DNAStringSet",
          dbFile = VignetteDB01,
          identifier = as.character(m1),
          verbose = TRUE)
}
names(genecalls) <- seq(length(genecalls))
x <- vapply(X = genecalls,
            FUN = function(x) {
              !is.null(x)
            },
            FUN.VALUE = vector(mode = "logical",
                               length = 1L),
            USE.NAMES = FALSE)
genecalls <- genecalls[x]

syn <- FindSynteny(dbFile = VignetteDB01,
                   verbose = TRUE)

save(syn,
     file = "~/Repos/SynExtend/data/syn.RData",
     compress = "xz")

save(genecalls,
     file = "~/Repos/SynExtend/data/genecalls.RData",
     compress = "xz")

###### -- NucleotideOverlap ---------------------------------------------------

linked_features <- NucleotideOverlap(SyntenyObject = syn,
                                     GeneCalls = genecalls,
                                     Verbose = TRUE)

save(linked_features,
     file = "~/Repos/SynExtend/data/linked_features.RData",
     compress = "xz")

###### -- PrepareSeqs ---------------------------------------------------------
# we're not using this function anymore, but we're creating a tmp db for the 
# rest of the example data because we want to ship the db without AAs appended
# for space considerations

system(command = paste("cp",
                       VignetteDB01,
                       VignetteDB02))

# PrepareSeqs(SynExtendObject = linked_features,
#             DataBase = VignetteDB02,
#             Verbose = TRUE)

###### -- PairSummaries -------------------------------------------------------

drv <- dbDriver("SQLite")
conn01 <- dbConnect(drv = drv,
                    VignetteDB02)

init_pairs <- SummarizePairs(SynExtendObject = linked_features,
                             DataBase = conn01,
                             Verbose = TRUE)

save(init_pairs,
     file = "~/Repos/SynExtend/data/init_pairs.RData",
     compress = "xz")

###### -- Evaluation ----------------------------------------------------------
# ensure this works, but no reason to save off the result in the example data

eval_pairs <- EvaluatePairs(InputPairs = init_pairs,
                            DataBase01 = conn01,
                            EvaluationMethod = "none",
                            FDRCriteria = c("Delta_Background" = 0.001))

###### -- Clustering ----------------------------------------------------------
# ensure this works, but no reason to save off the result in the example data

init_sets <- DisjointSet(Pairs = eval_pairs,
                         Verbose = TRUE)

# save(init_sets,
#      file = "~/Repos/SynExtend/data/init_sets.RData",
#      compress = "xz")

# Endosymbionts_Pairs02 <- ClusterByK(SynExtendObject = Endosymbionts_Pairs01,
#                                     ClusterScalar = 5,
#                                     ShowPlot = TRUE,
#                                     Verbose = TRUE)
# 
# save(Endosymbionts_Pairs02,
#      file = "~/Repos/SynExtend/data/Endosymbionts_Pairs02.RData",
#      compress = "xz")

###### -- BlockReconciliation -------------------------------------------------

# Endosymbionts_Pairs03 <- ExpandDiagonal(SynExtendObject = Endosymbionts_Pairs02[Endosymbionts_Pairs02$ClusterID %in% as.integer(names(which(attr(x = Endosymbionts_Pairs02,
#                                                                                                                                                  which = "Retain")))), ],
#                                         DataBase = CONN01,
#                                         Verbose = TRUE)
# save(Endosymbionts_Pairs03,
#      file = "~/Repos/SynExtend/data/Endosymbionts_Pairs03.RData",
#      compress = "xz")

###### -- DisjointSet ---------------------------------------------------------

# Endosymbionts_Sets <- DisjointSet(Pairs = Endosymbionts_Pairs03,
#                                   Verbose = TRUE)
# 
# save(Endosymbionts_Sets,
#      file = "~/Repos/SynExtend/data/Endosymbionts_Sets.RData",
#      compress = "xz")

###### -- ExtractBy ----------------------------------------------------------- 
# no functions in the pipeline beyond this function
# no need to save off this for examples
# Endosymbionts_Gene_Communities <- ExtractBy(x = Endosymbionts_Pairs03,
#                                             y = VignetteDB02,
#                                             z = Endosymbionts_Sets,
#                                             Verbose = TRUE)
# 
# save(Endosymbionts_Gene_Communities,
#      file = "~/Packages/SynExtend/data/Endosymbionts_Gene_Communities.RData",
#      compress = "xz")

###### -- CompetePairs --------------------------------------------------------

# no functions in the pipeline beyond this function
# so we don't need to save this for now...
# Endosymbionts_Pairs04 <- CompetePairs(SynExtendObject = Endosymbionts_Pairs01,
#                                       Verbose = TRUE)
# save(Endosymbionts_Pairs04,
#      file = "~/Packages/SynExtend/data/Endosymbionts_Pairs04.RData",
#      compress = "xz")



