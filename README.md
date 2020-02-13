# SynExtend

**Table of contents:**

*[Introduction](#introduction)
*[Installation](#installation)
*[Usage](#usage)
*[Reporting Bugs and Errors](#Reporting-Bugs-and-Errors)

## Introduction

SynExtend is a package of tools for working with synteny objects.

## Installation

SynExtend is currently in the submission process to Bioconductor and will be be available there upon acceptance, until then it can be installed from this repository via:

```r
library("devtools")
install_github(repo = "npcooley/SynExtend")
library("SynExtend")
```

Additionally, SynExtend is maintained in a Docker container at: Coming Soon!

## Usage

Currently SynExtend's major function is facilitating the prediction of orthologous gene pairs from synteny maps. The functions herein rely on synteny maps built by the `FindSynteny` function in the package `DECIPHER`. The process is relatively straightforward and only requires assemblies and genecalls, which can be conveniently pulled from the NCBI using entrez. However, as long as given assemblies have associated gene calls, or can have genecalls generated, they will work just fine.

Using entrez, this process would start like this:
```r
library(SynExtend)

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
                   timeout = 300L) # timeout argument is required for RStudio only

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
```

With these matched character vectors describing files on the NCBI FTP site, we can pull given assemblies, and their associated genecalls.

Assemblies are pulled using functions within the package `Biostrings` and placed in a sqlite database using the package `DECIPHER`. Note, `DBPATH` can be a specific local file location, though a tempfile works fine for this example.

```r
DBPATH <- tempfile()

for (m1 in seq_along(FNAs)) {
  X <- readDNAStringSet(FNAs[m1])
  X <- X[order(width(X),
               decreasing = TRUE)]
  Seqs2DB(seqs = X,
          type = "XStringSet",
          dbFile = DBPATH,
          identifier = as.character(m1),
          verbose = TRUE)
}
```

A function within the package `SynExtend` allows us to import GFF3 formatted gene annotations into a simple data.frame, and names are added to this generated list for later use.

```r
GeneCalls <- vector(mode = "list",
                    length = length(GFFs))

for (m1 in seq_along(GFFs)) {
  GeneCalls[[m1]] <- gffToDataFrame(GFF = GFFs[m1],
                                    Verbose = TRUE)
  
}

names(GeneCalls) <- seq(length(GeneCalls))
```

A synteny object is generated.

```r
Syn <- FindSynteny(dbFile = DBPATH)
```

Genecalls are overlaid on the synteny map.

```r
Links <- NucleotideOverlap(SyntenyObject = Syn,
                           GeneCalls = GeneCalls,
                           LimitIndex = FALSE,
                           Verbose = TRUE)
```

And positions on the grid of gene calls that are occupied by syntenic hits are recorded, and evaluated with a simple model to remove poor links. The function `PairSummaries` has arguments allowing a user align, and collect PIDs for each predicted pair, but that step is generally too cumbersome for a README example.

```r
PairedGenes <- PairSummaries(SyntenyLinks = Links,
                             GeneCalls = GeneCalls,
                             DBPATH = DBPATH,
                             PIDs = FALSE,
                             Model = "Global",
                             Verbose = TRUE)
```

A more self-contained example is maintained in this package's Bioconductor vignette, and can be found here: ...

## Reporting Bugs and Errors

Feel free to email me at npc19@pitt.edu, or submit bugs here on github.
