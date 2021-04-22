# SynExtend

**Table of contents:**
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Reporting Bugs and Errors](#Reporting-Bugs-and-Errors)

## Introduction

SynExtend is a package of tools for working with synteny objects generated from the Bioconductor Package DECIPHER.

## Installation

SynExtend is currently present in Bioconductor and the current release version can be installed via:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

BiocManager::install("SynExtend")
library("SynExtend")
```

The devel version in Bioconductor can be installed via:

```r
BiocManager::install("SynExtend", version = "devel")
library("SynExtend")
```

The version in this repository (which should always be up to date with devel - and sometimes ahead of release) can be installed with devtools via:

```r
devtools::install_github(repo = "npcooley/synextend")
library("SynExtend")
```

Additionally, SynExtend is maintained in a Docker container on [Dockerhub](https://hub.docker.com/repository/docker/npcooley/synextend) that is [here](https://github.com/npcooley/SynContainer). Only the current release version is maintained.

## Usage

Currently SynExtend's major function is facilitating the prediction of orthologous gene pairs from synteny maps. The functions herein rely on synteny maps built by the `FindSynteny()` function in the package `DECIPHER`. The process is relatively straightforward and only requires assemblies and genecalls, which can be conveniently pulled from the NCBI using [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/). However, as long as given assemblies have associated gene calls, or can have gene calls generated, they will work just fine. `DECIPHER` now contains the function `FindGenes()` that produces high quality de novo gene calls.

Using Entrez Direct, this process would start like this:
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

Assemblies are pulled using functions within the package `Biostrings` and placed in a sqlite database using the package `DECIPHER`. Note, `DBPATH` can be a specific local file location, though a tempfile works fine for this example. The tempfile in this example will be destroyed upon closure of the active R session, so specific file paths are preferable for most user activity.

```r
DBPATH <- tempfile()

for (m1 in seq_along(FNAs)) {
  X <- readDNAStringSet(FNAs[m1])
  # as of SynExtend 1.3.6 contig name matching implemented and the default.
  # Users can set AcceptContigNames = FALSE to prevent attempts at name matching,
  # but this will force SynExtend to assume that users have matched the ordering of contigs
  # in gene calls objects to their respective positions in the synteny objects.
  # X <- X[order(width(X),
  #              decreasing = TRUE)]
  Seqs2DB(seqs = X,
          type = "XStringSet",
          dbFile = DBPATH,
          identifier = as.character(m1),
          verbose = TRUE)
}
```

A few key quirks to this process are that genomes are given an integer identifier inside the `DECIPHER` database, and that same integer must be used as the name for the associated list position in the the object `GeneCalls` that is generated below. ~~Also SynExtend currently does not care about header names in FASTA files, meaning that the FASTA is ordered by sequence width before import to the data base.~~ **Contig names are now accepted as of SynExtend 1.3.6** Feature information in the GFF3 file is assigned an index and ordered by the width of its source sequence in the function `gffToDataFrame`.  This re-ordering may conflict with the natural order of contigs in fasta and therefore it is recommended to keep the default of `AcceptContigNames` as `TRUE`. Users also have the option of using the function `FindGenes` in `DECIPHER`, or using the `rtracklayer` function `import`.

```r
GeneCalls <- vector(mode = "list",
                    length = length(GFFs))

for (m1 in seq_along(GFFs)) {
  GeneCalls[[m1]] <- gffToDataFrame(GFF = GFFs[m1],
                                    Verbose = TRUE)
  
}

names(GeneCalls) <- seq(length(GeneCalls))
```

A synteny object is generated. See `?FindSynteny` in R for further reading.

```r
Syn <- FindSynteny(dbFile = DBPATH)
```

Genecalls are overlaid on the synteny map.

```r
Links <- NucleotideOverlap(SyntenyObject = Syn,
                           GeneCalls = GeneCalls,
                           LimitIndex = FALSE,
                           AcceptContigNames = TRUE,
                           Verbose = TRUE)
```

And positions on the grid of gene calls that are occupied by syntenic hits are recorded and evaluated with a simple model to identify unlikely pairs. The function `PairSummaries` has arguments allowing a user align and collect PIDs for each predicted pair, this step is cumbersome but `PIDs = TRUE` is currently the default.

```r
PairedGenes <- PairSummaries(SyntenyLinks = Links,
                             DBPATH = DBPATH,
                             PIDs = FALSE,
                             Model = "Generic",
                             AcceptContigNames =  TRUE,
                             Verbose = TRUE)
```

A more self-contained example is maintained in this package's Bioconductor vignette, and can be found [here](https://www.bioconductor.org/packages/release/bioc/html/SynExtend.html)

SynExtend currently has an extensive TODO list that includes:
1. Improving model performance for quickly identifying spuriously linked pairs.
2. Parsing communities of predicted pairs.
3. User friendly accessor functions to help with printing and plotting data.

If users have specific requests I would love to incorporate them.

## Reporting Bugs and Errors

Feel free to email me at npc19@pitt.edu or submit issues here on github.
