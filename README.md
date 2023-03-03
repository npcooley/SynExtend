# SynExtend

**Table of contents:**
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Reporting Bugs and Errors](#Reporting-Bugs-and-Errors)

## Introduction

SynExtend is a package of tools for working with synteny objects generated from the Bioconductor Package DECIPHER.

## Installation

SynExtend is present in the package repository Bioconductor and the current `release` version can be installed via:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

BiocManager::install("SynExtend")
library("SynExtend")
```

The `devel` version in Bioconductor can be installed via:

```r
BiocManager::install("SynExtend", version = "devel")
library("SynExtend")
```

The version in this repository (which should always be up to date with `devel` and frequently ahead of `release`) can be installed with devtools via:

```r
devtools::install_github(repo = "npcooley/synextend")
library("SynExtend")
```

Additionally, SynExtend is maintained in a Docker container on [Dockerhub](https://hub.docker.com/repository/docker/npcooley/synextend) that is [here](https://github.com/npcooley/SynContainer).

## Usage

### Prediction of Orthologous Gene Pairs

SynExtend facilitates the prediction of orthologous gene pairs from synteny maps. The functions herein rely on synteny maps built by the `FindSynteny()` function in the package `DECIPHER`. The process is relatively straightforward and only requires assemblies and genecalls, which can be conveniently pulled from the NCBI using [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/). However, as long as given assemblies have associated gene calls, or can have gene calls generated, they will work just fine. `DECIPHER` now contains the function `FindGenes()` that produces high quality de novo gene calls.

Using Entrez Direct, this process would start like this:
```r
library(SynExtend)

EntrezQuery <- paste("esearch -db assembly ",
                     "-query '",
                     "kitasatospora[organism] ",
                     'AND "complete genome"[filter] ', # only complete genomes
                     'AND "refseq has annotation"[properties] ', # only genomes with annotations
                     'AND "latest refseq"[filter] ', # only latest
                     "NOT anomalous[filter]' ",
                     '| ',
                     'esummary ',
                     '| ',
                     'xtract -pattern DocumentSummary -element FtpPath_RefSeq',
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

With these matched character vectors describing files on the NCBI FTP site, we can pull given assemblies, and their associated `GFF` files.

Assemblies are pulled using functions within the package `Biostrings` and placed in a sqlite database using the package `DECIPHER`. Note, `DBPATH` can be a specific local file location, though a tempfile works fine for this example. The tempfile in this example will be destroyed upon closure of the active R session, so specific file paths are preferable for most user activity. The code chunk below provides an example of pulling down assemblies from the NCBI and using `DECIPHER`'s built in gene finder to generate gene calls.

```r
DBPATH <- tempfile()

GC01 <- vector(mode = "list",
               length = length(FNAs))

for (m1 in seq_along(FNAs)) {
  X <- readDNAStringSet(FNAs[m1])
  Seqs2DB(seqs = X,
          type = "XStringSet",
          dbFile = DBPATH,
          identifier = as.character(m1),
          verbose = TRUE)
  
  GC01[[m1]] <- FindGenes(myDNAStringSet = X)
}

names(GC01) <- seq(length(FNAs))
```

A current quirk to this process is that genecalls must be present in a named list, where the names can be matched to integer identifiers in the `DECIPHER` database used to construct a synteny object. Users have the option of using the function `FindGenes` in `DECIPHER`, using the `rtracklayer` function `import` to import a valid `GFF` file, or importing a `GFF` with the `SynExtend` function `gffToDataFrame` which performs some data munging to convert the `GFF` into a `DataFrame` for user convenience.

```r
GC02 <- vector(mode = "list",
               length = length(GFFs))

for (m1 in seq_along(GFFs)) {
  GC02[[m1]] <- gffToDataFrame(GFF = GFFs[m1],
                               Verbose = TRUE)
  
}

names(GC02) <- seq(length(GC02))
```

With assemblies and genecalls now in a user's workspace, users can build a synteny map with `FindSynteny`. See `?FindSynteny` in R for further reading. Finding and tabulating where genes are connected by syntenic k-mers is then performed with `NucleotideOverlap`, and that raw data is converted into a convenient `data.frame` with the function `PairSummariers`.

```r
Syn <- FindSynteny(dbFile = DBPATH)

Links <- NucleotideOverlap(SyntenyObject = Syn,
                           GeneCalls = GC01,
                           Verbose = TRUE)
                           
Pairs01 <- PairSummaries(SyntenyLinks = Links,
                         PIDs = TRUE,
                         Verbose = TRUE)
                             
Clusters01 <- DisjointSet(Pairs = Pairs01,
                          Verbose = TRUE)
```

Pairs where predictions create conflicts can be removed by a simple subsetting function.

```r
Pairs02 <- SubSetPairs(CurrentPairs = Pairs01,
                       Verbose = TRUE)

Clusters02 <- DisjointSet(Pairs = Pairs02,
                          Verbose = TRUE)
```


A more self-contained example is maintained in this package's Bioconductor vignette, and can be found [here](https://www.bioconductor.org/packages/release/bioc/html/SynExtend.html)

### Creating Functional Association Networks

SynExtend also allows for prediction of functional association networks from evidence of correlated selective pressures between of pairs of clusters of orthologous genes (COGs). Given a list of gene trees or presence/absence data, the `ProtWeaver` class in SynExtend uses multiple algorithms to compute likelihood of pairwise functional association. This list can be generated with `DisjointSet(...)`, or can be created by the user. 

An example of using `ProtWeaver` can be done using the built-in example *Streptomyces* dataset. This will predict functional association of 200 genes from 301 genomes. 

```r
exData <- get(data("ExampleStreptomycesData"))
pw <- ProtWeaver(exData$Genes)

# For faster runtime, set subset as follows:
# predictions <- predict(pw, mySpeciesTree=exData$Tree, subset=1:20)
predictions <- predict(pw, mySpeciesTree=exData$Tree)

# View number of genes and number of predictions:

# Print out the adjacency matrix:
predictions

# Print out pairwise associations:
as.data.frame(predictions)

# Plot genes using force-directed embedding:
plot(predictions)
```

Current algorithms used for prediction are Jaccard distance, Hamming distance, Mutual Information, Direct Coupling Analysis, Phylogenetic Gain/Loss events, Co-localization, MirrorTree, and ContextTree. See `?ProtWeaver` for more details.

### Planned Updates

SynExtend currently has an extensive TODO list that includes:
1. Parsing communities of predicted pairs.
2. User friendly accessor functions to help with printing and plotting data.

If users have specific requests I would love to incorporate them.

## Reporting Bugs and Errors

Feel free to email me at npc19@pitt.edu or submit issues here on GitHub. Questions related to `ProtWeaver` can be sent to ahl27@pitt.edu or submitted similarly on GitHub.
