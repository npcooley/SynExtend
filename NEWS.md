# SynExtend 1.3.8

* `PairSummaries` now allows users to fill in specific matching gaps in blocks of predicted pairs with the arguments `AllowGaps` and `OffSetsAllowed`.

# SynExtend 1.3.7

* Adjustments to progress bars in both `PairSummaries` and `NucleotideOverlap`.
* PID prediction models in `PairSummaries` adjusted.

# SynExtend 1.3.6

* Contig name matching has been implemented. Scheme expects users to follow NCBI contig naming and gff formats, accepting contig names from gffs directly, and removing the first whitespace and everything thereafter from FASTA headers. Contig name matching can be disabled if users wish, using the argument `AcceptContigNames`, but ensuring that the correct contigs in GeneCalls objects are matched to the appropriate contigs in Synteny objects are then the user's responsibility.

# SynExtend 1.3.5

* `PairSummaries` now translates sequences based on `transl_table` attributes provided by gene calls
* `PairSummaries` now uses a generic model for predicting PID
* `gffToDataFrame` now parses out the `transl_table` attribute

# SynExtend 1.3.2

* Minor changes to `NucleotideOverlap`
* Major changes to `PairSummaries` - can now take in objects of class `Genes` build by the DECIPHER function `FindGenes()`

# SynExtend 0.99.30

* Vignette and help files edited for clarity

# SynExtend 0.99.1

* `SynExtend` submitted to Bioconductor
* Added function `gffToDataFrame`
* Added function `NucleotideOverlap`
* Added function `PairSummaries`